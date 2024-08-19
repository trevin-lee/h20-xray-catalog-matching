from pathlib import Path
import yaml

class ConfigLoader:

    def __init__(self, config_path: str):
        config_path: Path = Path(config_path).resolve()
        try:
            config_file = open(str(config_path), 'r')
            self.yaml_dict = yaml.safe_load(config_file)
            self.load_member_vars()
        except FileNotFoundError:
            print(f"{config_file} not found")
        except yaml.YAMLError:
            print(f"Error loading YAML file: {yaml.YAMLError}")


    def load_member_vars(self):
        data = self.yaml_dict
        
        self.EROSITA_PATH = data['eROSITA_file_path']
        self.EDFF_PATH = data['EDFF_file_path']
        self.EDFN_PATH = data['EDFN_file_path']
        self.CDFS_PATH = data['CDFS_file_path']

        self.EDFF_CDFS_PATH = data['edff_cdfs_path']
        self.EDFF_EROSITA_PATH = data['edff_erosita_path']

        self.EDFF_SED_DATA = data['edff_sed_data']

        self.SAVE_DATA_PATH = data['save_data_path']
        self.RAW_DATA_PATH = data['raw_data_path']
        self.INPUT_CATALOG_PATH = data['input_catalogs_path']

        self.ADJ_FLUX_RAD = data['adj_flux_rad']
        self.MATCH_SEARCH_RAD = data['match_search_rad']
        self.BACKGROUND_INNER_RAD = data['background_inner_rad']
        self.BACKGROUND_OUTER_RAD = data['background_outer_rad']

        self.INTERPOLATION_BINS = data['interpolation_bins']