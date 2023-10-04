import sys
import os
import enum
import math
import numpy as np
import typing as ty

# pj imports
import phylojunction.utility.helper_functions as pjh

__author__ = "Fabio K. Mendes"
__email__ = "f.mendes@wustl.edu"


class GeoFeatureRelationship(enum.Enum):
    WITHIN = "within"
    BETWEEN = "between"


class GeoFeatureType(enum.Enum):
    CATEGORICAL = "categorical"
    QUANTITATIVE = "quantitative"


class GeoFeature():
    time_idx: int  # epoch index
    feat_idx: int  # should be unique to this feature
    n_regions: int  # automatically computed
    feat_rel: GeoFeatureRelationship
    feat_type: GeoFeatureType
    value: np.array
    name: str

    def __init__(self,
                 time_idx: int,
                 feat_idx: int,
                 feat_rel: GeoFeatureRelationship,
                 feat_type: GeoFeatureType,
                 feat_vals: ty.Optional[np.array] = None,
                 feat_name: ty.Optional[str] = ""):

        self.time_idx = time_idx
        self.feat_idx = feat_idx
        self.feat_rel = feat_rel
        self.feat_type = feat_type
        self.value = feat_vals
        self.name = feat_name
        self.n_regions = None

        if feat_type == GeoFeatureRelationship.BETWEEN:
            self.ndim = 2
        
        if self.value is not None:
            self._initialize_n_regions()

        self._initialize_str_representation()

    # init methods
    def _initialize_str_representation(self) -> None:
        self.str_representation = \
            "Feature | " + self.feat_rel.value + "-" \
            + self.feat_type.value + " " + str(self.feat_idx) + " | " \
            + "Epoch " + str(self.time_idx) 
        
        if self.name:
            str_replacement = "Feature (" + self.name + ")"
            self.str_representation = \
                self.str_representation.replace("Feature", str_replacement)

        if self.n_regions:
            self.str_representation += " | " + str(self.n_regions) \
                + " regions"

        if self.value is not None:
            if self.feat_rel == GeoFeatureRelationship.WITHIN:
                self.str_representation += "\n" + np.array2string(self.value)

            elif self.feat_rel == GeoFeatureRelationship.BETWEEN:
                self.str_representation += "\n" + np.array2string(
                    self.value.reshape(self.n_regions,
                                       self.n_regions))

    # setters
    def _set_attr_values(self, attr_name: str, val:
                         ty.Tuple[int,
                                  GeoFeatureRelationship,
                                  GeoFeatureType,
                                  str,
                                  np.array]) -> None:
        self.__setattr__(attr_name, val)

        # automatically setting self.n_regions
        if attr_name == "value":
            self._initialize_n_regions()

        # force update
        self._initialize_str_representation()

    # getters
    def _get_val_from_to(from_idx: int, to_idx: int) -> float:
        pass

    # health
    def _initialize_n_regions(self) -> None:
        if self.feat_rel == GeoFeatureRelationship.WITHIN:
            self.n_regions = len(self.value)

        elif self.feat_rel == GeoFeatureRelationship.BETWEEN:
            self.n_regions = int(math.sqrt(len(self.value)))
    
    # other
    def __str__(self) -> None:
        return self.str_representation


class FeatureCollection():

    geofeat_list: ty.List[GeoFeature] = list()

    # ways in which to store (for later grabbing) features
    #
    # (1) by name (assumes name cannot be repeated across feature
    # types and/or relationships)
    #
    # TODO: add check here!!!!
    feat_name_epochs_dict: ty.Dict[str, ty.Dict[int, GeoFeature]] \
        = pjh.autovivify(2)
    # (2) by feature type, relationship and feature number (index)
    feat_type_rel_featid_epochs_dict: \
        ty.Dict[GeoFeatureType, 
                ty.Dict[GeoFeatureRelationship,
                        ty.Dict[int,
                                ty.Dict[int, GeoFeature]
                                ]]] \
                                = pjh.autovivify(4)

    region_name_idx_dict: ty.Dict[str, int] = dict()
    region_idx_name_dict: ty.Dict[int, str] = dict()

    def __init__(self, feat_summary_fp: str) -> None:
        self._check_filepaths(feat_summary_fp)

        # initializes:
        #     region_name_idx_dict
        #     region_idx_name_dict
        self._read_feat_summary_init_feats(feat_summary_fp)

        # initializes self.feat_name_epochs_dict
        self._init_feat_name_epochs_dict()

        # initializes self.feat_type_rel_featid_epochs_dict
        self._init_feat_type_rel_featid_epochs_dict()

        # debugging
        # self._debug_print_dicts()

    # init methods
    def _read_feat_summary_init_feats(
            self,
            feat_summary_fp):
        """Populate feature list with initial info
        
        Open, read and parse feature summary .csv file, and populate
        class member list of features with initial information
        """

        with open(feat_summary_fp, "r") as infile:
            header = infile.readline()
            if header != ("time_index,feature_index,feature_relationship,"
                          "feature_type,feature_path,name\n"):
                exit(("ERROR: Expecting header:\ntime_index,feature_index,"
                      "feature_relationship,feature_type,feature_path,name\n"))
            
            # linecount = 0
            # feat_unique_names = set([])
            for line in infile:
                line = line.rstrip()
                time_idx, feat_idx, feat_rel, feat_type, feat_fp, feat_name \
                    = line.split(",")

                # feat_unique_names.add(feat_name)

                geofeat = GeoFeature(
                    int(time_idx),
                    int(feat_idx),
                    (GeoFeatureRelationship.WITHIN \
                    if feat_rel == "within" \
                    else GeoFeatureRelationship.BETWEEN),
                    (GeoFeatureType.CATEGORICAL \
                    if feat_rel == "categorical" \
                    else GeoFeatureType.QUANTITATIVE),
                    feat_name = feat_name)
                
                self._check_filepaths(feat_fp)
                
                with open(feat_fp, "r") as infile2:
                    regions_str = infile2.readline().rstrip()
                    regions_list = regions_str.split(",")

                    for idx, reg in enumerate(regions_list):
                        self.region_name_idx_dict[reg] = idx
                        self.region_idx_name_dict[idx] = reg

                    vals_list = list()
                    for line in infile2:
                        line = line.rstrip().split(",")
                        vals_list.extend(float(v) for v in line)
                    
                    vals_array = np.array(vals_list)
                
                geofeat._set_attr_values("value", vals_array)
                
                # debugging
                # print(geofeat)

                self.geofeat_list.append(geofeat)

                # linecount += 1
            
            # if len(feat_unique_names) != linecount:
            #     exit("ERROR: Detected " + str(len(feat_unique_names)) \
            #          + " unique feature names, but summary file had " \
            #          + str(linecount) + " lines minus the header. Exiting.")

    def _init_feat_name_epochs_dict(self):
        for geofeat in self.geofeat_list:
            self.feat_name_epochs_dict[geofeat.name][geofeat.time_idx] \
                = geofeat

    def _init_feat_type_rel_featid_epochs_dict(self):
        for geofeat in self.geofeat_list:
            fty = geofeat.feat_type
            fr = geofeat.feat_rel
            fid = geofeat.feat_idx
            fti = geofeat.time_idx
            self.feat_type_rel_featid_epochs_dict[fty][fr][fid][fti] \
                = geofeat

    # health check methods
    def _check_filepaths(self,
                         feat_summary_fp) -> None:
        """Check files exist in initialization file paths"""

        if not os.path.isfile(feat_summary_fp):
            exit("Could not find " + feat_summary_fp + ". Exiting.")

    # getters
    def get_feat_by_name(self, feat_name: str, time_idx: ty.Optional[int] = None) \
            -> ty.Tuple[float, ty.List[float]]:

        feat_dict_all_epochs = self.feat_name_epochs_dict[feat_name]

        if time_idx:
            return feat_dict_all_epochs[time_idx]

        else:
            return [v for k, v in feat_dict_all_epochs.items()]

    # debug
    def _debug_print_dicts(self) -> None:
        for feat_name, feat_dict in self.feat_name_epochs_dict.items():
            for t, feat in feat_dict.items():
                print(feat_name, ":\n", feat, sep="")
        
        for fty, feat_dict1 in self.feat_type_rel_featid_epochs_dict.items():
            for fr, feat_dict2 in feat_dict1.items():
                for fid, feat_dict3 in feat_dict2.items():
                    for fti, feat in feat_dict3.items():
                        print(fty.value + "-" + fr.value + " feature " \
                              + str(fid) + " epoch " + str(fti),
                              feat,
                              sep="\n")


# callable type for type hinting requirement functions
class MyCallableType(ty.Protocol):
    def __call__(self, feat_col: FeatureCollection) -> None:
        ...

class FeatureQuery():
    feat_coll: FeatureCollection

    def __init__(self, feat_coll) -> None:
        self.feat_coll = feat_coll

    # requirement function 1
    @classmethod
    def cw_features_are_one(cls,
                            feat_col: FeatureCollection,
                            feat_name: ty.List[str] = [],
                            feat_id: ty.List[int] = []) -> ty.List[float]:
    
        if feat_name:
            pass

        elif feat_id:
            pass

        else:
            exit(("ERROR: Function \'cw_features_are_one\' needs either one "
                    "or more categorical-within feature names or feature ids."
                    ". Exiting."))
            
        return [1.0]
    
    # # requirement function 2
    @classmethod
    def qw_features_threshold(cls,
                            feat_col: FeatureCollection,
                            threshold: float,
                            feat_name: ty.List[str] = [],
                            feat_id: ty.List[int] = []) -> ty.List[float]:
    
        if feat_name:
            pass

        elif feat_id:
            pass

        else:
            exit(("ERROR: Function \'qw_features_threshold\' needs either one "
                    "or more categorical-within feature names or feature ids."
                    ". Exiting."))
            
        return [2.0]
        
    def when_has_it_happened(self,
                             geo_event_name: str,
                             requirement_fn: MyCallableType) \
            -> ty.List[float]:
        """Return list of times when events happened
        
        Args:
            geo_event_name (str): the name one is giving the event
            requirement_fn (function): this function receives a feature
                collection and a list of integers defining the involved
                regions, and returns a list of times

        Returns:
            list (float): time(s) the event took place
        """

        return requirement_fn


if __name__ == "__main__":
    # sys.argv[1]: feature summary file path

    fc = FeatureCollection(sys.argv[1])

    # testing
    # print(fc.get_feat_by_name("qb_1", time_idx=1))
    # feat_list = fc.get_feat_by_name("qb_1")
    # for feat in feat_list:
    #     print(feat)

    fq = FeatureQuery(fc)
    
    requirement_fn1 = \
        FeatureQuery.cw_features_are_one(fc, "cw1")
    
    res1 = fq.when_has_it_happened("barrier", requirement_fn1)
    print(res1)

    thresh = 0.5
    requirement_fn2 = \
        FeatureQuery.qw_features_threshold(fc, "qw1", thresh)
    
    res2 = fq.when_has_it_happened("barrier", requirement_fn2)
    print(res2)

    # let's say we want some weird function that is not
    # defined inside FeatureQuery...
    #
    # we want, say, that a categorical between feature == 1
    # and then that a quantitative between feature is < threshold
    def custom_requirement_fn(feature_collection,
                              cb_feat_name,
                              qb_feat_name,
                              threshold):
        
        if not (cb_feat_name and qb_feat_name):
            exit(("ERROR: Function \'qw_features_threshold\' needs either one "
                    "or more categorical-within feature names or feature ids."
                    ". Exiting."))

        else:
            pass
            
        return [3.0]

    requirement_fn3 = custom_requirement_fn(fc,
                                            "cb1",
                                            "qb1",
                                            thresh)
    res3 = fq.when_has_it_happened("barrier", requirement_fn3)
    print(res3)