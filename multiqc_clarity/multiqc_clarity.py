import csv
import traceback
from collections import OrderedDict, defaultdict
from os.path import isfile, join

from multiqc.utils import report, config
from multiqc.modules.base_module import BaseMultiqcModule
from multiqc.plots import table

import logging


class MultiQC_clarity_metadata(BaseMultiqcModule):
    def __init__(self):

        self.log = logging.getLogger('multiqc')

        # Check that this plugin hasn't been disabled
        if config.kwargs.get('disable_clarity', False) is True:
            self.log.info("Skipping MultiQC_Clarity as disabled on command line")
            return
        if getattr(config, 'disable_clarity', False) is True:
            self.log.debug("Skipping MultiQC_Clarity as specified in config file")
            return

        super(MultiQC_clarity_metadata, self).__init__(name='Clarity LIMS', anchor='clarity')

        self.intro = '''<p>The <a href="https://github.com/MultiQC/MultiQC_Clarity" target="_blank">MultiQC_Clarity</a>
            plugin fetches data from a specified
            <a href="https://www.genologics.com/clarity-lims/" target="_blank">Basespace Clarity LIMS</a> instance.</p>'''

        try:
            from genologics.lims import Lims
            from genologics import config as genologics_config
        except:
            self.log.warning("Importing genologics failed: " + traceback.format_exc())
            return
        
        try:
            BASEURI, USERNAME, PASSWORD, VERSION, MAIN_LOG = genologics_config.load_config(specified_config=config.kwargs.get('clarity_config'))
        except SystemExit:
            self.log.warning("Genologics config file is not specified as --clarity_config or in ~/.genologicsrc. "
                             "Skip running Clarity module")
            return
        
        self.lims = Lims(BASEURI, USERNAME, PASSWORD)
        self.metadata = {}
        self.header_metadata = {}
        self.general_metadata = {}
        self.tab_metadata = {}
        self.samples = []

        self.schema = getattr(config, 'clarity', None)
        if self.schema is None:
            self.log.debug("No config found for MultiQC_Clarity")
            return

        self.get_samples()
        if 'report_header_info' in self.schema:
            self.get_metadata('report_header_info')
        if 'general_stats' in self.schema:
            self.get_metadata('general_stats')
        if 'clarity_module' in self.schema:
            self.get_metadata('clarity_module')
        self.update_multiqc_report()
        self.make_sections()
        report.modules_output.append(self)

    def csv_file_from_samplesheet(self, sample_sheet):
        csv_lines = []
        with open(sample_sheet) as f:
            found_data = False
            for line in f:
                if found_data:
                    csv_lines.append(line.strip())
                else:
                    if line.strip().startswith('[Data]'):
                        found_data = True
        return csv_lines

    def get_raw_sample_names(self, csv_fpath, names):
        raw_sample_names = dict()
        with open(csv_fpath) as f:
            csv_reader = csv.DictReader(f)
            name_col = csv_reader.fieldnames[0]
            for r in csv_reader:
                correct_name = r['description'] if 'description' in r else r[name_col]
                if correct_name not in names:
                    continue
                raw_sample_names[correct_name] = r[name_col]
        return raw_sample_names

    def correct_sample_name(self, name):
        import re
        name = re.sub(r'_S\d+$', '', name)
        return name.replace('-', '_').replace('.', '_')

    def search_by_samplesheet(self, names):
        sample_sheet_fpath = config.kwargs['samplesheet']
        samples_by_container = defaultdict(dict)
        raw_names = dict((name, name) for name in names)
        if config.kwargs.get('bcbio_csv') and isfile(config.kwargs.get('bcbio_csv')):
            raw_names = self.get_raw_sample_names(config.kwargs['bcbio_csv'], names)

        correct_sample_names = dict((self.correct_sample_name(raw_names[name]), name) for name in names)
        for row in csv.DictReader(self.csv_file_from_samplesheet(sample_sheet_fpath), delimiter=','):
            sample_name = row['SampleName']
            if sample_name not in correct_sample_names.keys():
                continue
            sample_id = row['SampleID'] if 'SampleID' in row else row['Sample_ID']
            container, sample_well = row['SamplePlate'], row['SampleWell'].replace('_', ':')
            sample_artifacts = self.lims.get_artifacts(samplelimsid=sample_id)
            if sample_artifacts:
                sample = sample_artifacts[0].samples[0]
                sample.name = correct_sample_names[sample_name]
                self.samples.append(sample)
            else:
                samples_by_container[container][sample_well] = sample_name

        for container_id, samples in samples_by_container.items():
            artifacts = self.lims.get_artifacts(containerlimsid=container_id)
            if not artifacts:
                continue
            placements = artifacts[0].container.get_placements()
            for well, sample_name in samples.items():
                sample = placements[well].samples[0]
                sample.name = correct_sample_names[sample_name]
                self.samples.append(sample)

    def get_samples(self):
        if config.kwargs.get('clarity_project_name'):
            pj = self.lims.get_projects(name=config.kwargs['clarity_project_name'])
            self.samples = pj.samples
            self.log.info("Found {} in LIMS.".format(config.kwargs['clarity_project_name']))
        else:
            names = set()
            for x in report.general_stats_data:
                names.update(x.keys())
            for d in report.saved_raw_data.values():
                try:
                    self.names.update(d.keys())
                except AttributeError:
                    pass
            if not config.kwargs.get('clarity_skip_edit_names'):
                names = self.edit_names(names)

            self.log.debug("Looking into Clarity for samples {}".format(", ".join(names)))
            if config.kwargs.get('samplesheet'):
                self.search_by_samplesheet(names)
            if not self.samples:
                try:
                    for name in names:
                        matching_samples = self.lims.get_samples(name=name)
                        if not matching_samples:
                            self.log.error("Could not find a sample matching {0}, skipping.".format(name))
                            continue
                        if len(matching_samples) > 1:
                            self.log.error("Found multiple samples matching {0}, skipping".format(name))
                            continue
                        self.samples.append(matching_samples[0])
                except Exception as e:
                    self.log.warn("Could not connect to Clarity LIMS: {}".format(e))
                    return None
            self.log.info("Found {} out of {} samples in LIMS.".format(len(self.samples), len(names)))


    def edit_names(self, names):
        edited=[]
        for name in names:
            if name.endswith("_1") or name.endswith("_2"):
                edited.append(name[:-2])
            elif name.endswith("_R1") or name.endswith("_R2"):
                edited.append(name[:-3])
            else:
                edited.append(name)

        return edited

    def flatten_metadata(self, metadata):
        for first_level in metadata:
            for second_level in metadata[first_level]:
                if isinstance(metadata[first_level][second_level], set) or isinstance(metadata[first_level][second_level], list):
                    metadata[first_level][second_level] = ", ".join(metadata[first_level][second_level])

        return metadata

    def get_project_metadata(self, udfs):
        project_metadata={}
        for sample in self.samples:
            project_metadata[sample.project.name]={}
            for udf in udfs:
                if udf in sample.project.udf:
                    try:
                        project_metadata[sample.project.name][udf].add(str(sample.project.udf[udf]))
                    except:
                        project_metadata[sample.project.name][udf] = set()
                        project_metadata[sample.project.name][udf].add(str(sample.project.udf[udf]))

        return self.flatten_metadata(project_metadata)

    def get_sample_metadata(self, udfs):
        sample_metadata={}
        for sample in self.samples:
            sample_metadata[sample.name]={}
            for udf in udfs:
                if udf in sample.udf:
                    try:
                        sample_metadata[sample.name][udf].add(str(sample.udf[udf]))
                    except:
                        sample_metadata[sample.name][udf] = set()
                        sample_metadata[sample.name][udf].add(str(sample.udf[udf]))
            sample_type = sample_metadata[sample.name].pop('Sample Tissue') if \
                'Sample Tissue' in sample_metadata[sample.name] else sample_metadata[sample.name].pop('Sample Type')
            sample_link = join(self.lims.baseuri, 'clarity', 'search?scope=Sample&query=' + sample.id)
            sample_metadata[sample.name]['Sample Type'] = '<a href="' + sample_link + '" target="_blank">' + sample_type.pop() + '</a>'
            report.lims_added = True
        return self.flatten_metadata(sample_metadata)


    def get_metadata(self, part):
        for key in self.schema[part]:
            if key == 'Project':
                metadata = self.get_project_metadata(self.schema[part]['Project'])
            elif key == 'Sample':
                metadata = self.get_sample_metadata(self.schema[part]['Sample'])
            else:
                metadata = self.get_artifact_metadata(self.schema[part])

            if part == "report_header_info":
                self.header_metadata.update(metadata)
            elif part == "general_stats":
                self.general_metadata.update(metadata)
            else:
                self.tab_metadata.update(metadata)


    def get_artifact_metadata(self, pt_to_udfs):
        artifact_metadata={}
        for sample in self.samples:
            artifact_metadata[sample.name]={}
            for process_type in pt_to_udfs:
                if process_type == 'Sample':
                    continue
                if process_type == 'Project':
                    continue
                artifacts = self.lims.get_artifacts(sample_name=sample.name, process_type=process_type)
                for udf_name in pt_to_udfs[process_type].get("outputs", []):
                    values = []
                    for artifact in artifacts:
                        if udf_name in artifact.udf:
                            values.append(str(artifact.udf[udf_name]))

                    artifact_metadata[sample.name][udf_name]=values

                processes = set([art.parent_process for art in artifacts])
                inputs=[]
                for p in processes:
                    inputs.extend([art for art in p.all_inputs() if sample.name in [s.name for s in art.samples]])
                for udf_name in pt_to_udfs[process_type].get("inputs", []):
                    values = []
                    for artifact in inputs:
                        if udf_name in artifact.udf:
                            values.append(str(artifact.udf[udf_name]))

                    artifact_metadata[sample.name][udf_name]=values

        return self.flatten_metadata(artifact_metadata)


    def update_multiqc_report(self):
        if config.report_header_info is None:
            config.report_header_info = []
        for first_level in self.header_metadata:
            d = {}
            for key in self.header_metadata[first_level]:
                d[key] = self.header_metadata[first_level][key]
            config.report_header_info.append(d)

        headers = {}
        for first_level in self.schema["general_stats"]:
            for header in self.schema["general_stats"][first_level]:
                headers[header] = {}
                if isinstance(self.schema["general_stats"][first_level][header], dict):
                    for subsubkey, cfg in self.schema["general_stats"][first_level][header].items():
                        if subsubkey == 'multiply_by':
                            mby = str(cfg)[:]
                            headers[header]['modify'] = lambda x: float(x) * float(mby)
                        else:
                            headers[header][subsubkey] = cfg
                headers[header]['description'] = headers[header].get('description', '{} - {}'.format(first_level, header))
                headers[header]['namespace'] = headers[header].get('namespace', 'Clarity LIMS')
                headers[header]['scale'] = headers[header].get('scale', 'YlGn')

        report.general_stats_headers.append(headers)
        report.general_stats_data.append(self.general_metadata)

    def make_sections(self):
        headers = OrderedDict()
        for first_level in self.tab_metadata:
            for header in self.tab_metadata[first_level]:
                desc = header
                if header not in headers:
                    headers[header] = {}
                    for key in self.schema['clarity_module']:
                        if header in self.schema['clarity_module'][key]:
                            desc = key
                        elif isinstance(self.schema['clarity_module'][key], dict):
                            for subkey, val in self.schema['clarity_module'][key].items():
                                # print(val)
                                if val is None:
                                    break
                                elif header in val:
                                    desc = key
                                    if isinstance(val[header], dict):
                                        for subsubkey, cfg in val[header].items():
                                            if subsubkey == 'multiply_by':
                                                mby = str(cfg)[:]
                                                headers[header]['modify'] = lambda x: float(x) * float(mby)
                                            else:
                                                headers[header][subsubkey] = cfg

                    headers[header]['namespace'] = headers[header].get('namespace', desc)
                    headers[header]['title'] = headers[header].get('title', header)
                    headers[header]['description'] = headers[header].get('description', header)

        self.intro += table.plot(self.tab_metadata, headers)


