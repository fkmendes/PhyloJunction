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
            self.n_dim = 2
        
        if self.value is not None:
            # initialized here or when self.value is 
            # set by setter
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
    def get_val_from_to(self, from_idx: int, to_idx: int) -> float:
        if not self.n_regions:
            exit("ERROR: Feature " + self.name + " is not aware of " \
                 + "how many regions the feature is scored for. Exiting.")
        
        return self.value[from_idx * self.n_regions + to_idx]

    # health
    def _initialize_n_regions(self) -> None:
        if self.feat_rel == GeoFeatureRelationship.WITHIN:
            self.n_regions = len(self.value)

        elif self.feat_rel == GeoFeatureRelationship.BETWEEN:
            self.n_regions = int(math.sqrt(len(self.value)))
    
    # other
    def __str__(self) -> None:
        return self.str_representation


class GeoFeatureCollection():

    geofeat_list: ty.List[GeoFeature] = list()
    epoch_age_end_list_young2old: ty.List[float] = [0.0]  # present epoch ends at 0.0
    epoch_mid_age_list_young2old: ty.List[float] = [-math.inf]  # oldest epoch midpoint 
    is_timehet: bool = False

    # ways in which to store (for later grabbing) features
    #
    # (1) by name (assumes name cannot be repeated across feature
    # types and/or relationships); final dict is k: time_idx, v: feat
    #
    # TODO: add check here!!!!
    feat_name_epochs_dict: ty.Dict[str, ty.Dict[int, GeoFeature]] \
        = pjh.autovivify(2)
    # (2) by feature type, relationship and feature number (index);
    # final dict is k: time_idx, v: feat
    feat_type_rel_featid_epochs_dict: \
        ty.Dict[GeoFeatureType, 
                ty.Dict[GeoFeatureRelationship,
                        ty.Dict[int,
                                ty.Dict[int, GeoFeature]
                                ]]] \
                                = pjh.autovivify(4)

    region_name_idx_dict: ty.Dict[str, int] = dict()
    region_idx_name_dict: ty.Dict[int, str] = dict()

    def __init__(self,
                 feat_summary_fp: str,
                 age_summary_fp: str = "") -> None:
        
        self._check_filepaths(feat_summary_fp)

        # if only one epoch, age_summary not necessary
        if age_summary_fp:
            self._check_filepaths(age_summary_fp)

        # initializes:
        #     is_timehet
        #     region_name_idx_dict
        #     region_idx_name_dict
        self._read_feat_summary_init_feats(feat_summary_fp)

        # initializes self.feat_name_epochs_dict
        self._init_feat_name_epochs_dict()

        # initializes self.feat_type_rel_featid_epochs_dict
        self._init_feat_type_rel_featid_epochs_dict()

        # debugging
        # self._debug_print_dicts()

        # initializes
        #     self.epoch_age_end_list
        #     self.epoch_mid_age_list
        self._read_age_summary(age_summary_fp)

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

                if int(time_idx) > 1:
                    self.is_timehet = True

                geofeat = GeoFeature(
                    int(time_idx),
                    int(feat_idx),
                    (GeoFeatureRelationship.WITHIN \
                    if feat_rel == GeoFeatureRelationship.WITHIN.value \
                    else GeoFeatureRelationship.BETWEEN),
                    (GeoFeatureType.CATEGORICAL \
                    if feat_type == GeoFeatureType.CATEGORICAL.value \
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
                        vals_list.extend(int(v) \
                                         if feat_type == GeoFeatureType.CATEGORICAL.value \
                                         else float(v) \
                                         for v in line)
                    
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

    def _read_age_summary(
            self,
            age_summary_fp):
        """Populate epoch age end list
        
        Open, read and parse age summary .csv file, and populate
        class member list of epoch age ends
        """

        with open(age_summary_fp, "r") as infile:
            header = infile.readline()
            
            if header != ("index,age_end\n"):
                exit(("ERROR: Expecting header:\nindex,age_end\n"))

            content = infile.readlines()
            n_epochs_minus_one = len(content)
            epoch_age_ends_from_file = \
                [-1.0 for i in range(n_epochs_minus_one)]

            for line in content:
                idx, age_end = line.rstrip().split(",")
                # note the offset!
                epoch_age_ends_from_file[int(idx)-1] = float(age_end)

        # actual initialization
        self.epoch_age_end_list_young2old += epoch_age_ends_from_file

        # debugging
        # print(self.epoch_age_end_list_young2old)

        n_epochs = len(self.epoch_age_end_list_young2old)
        mid_ages_minus_last_epoch = list()
        for idx in range(len(self.epoch_age_end_list_young2old)):
            if idx < n_epochs-1:
                younger_age = self.epoch_age_end_list_young2old[idx]
                older_age = self.epoch_age_end_list_young2old[idx+1]
                mid_age = (younger_age + older_age) / 2.0

                mid_ages_minus_last_epoch.append(mid_age)

                # debugging
                # print(younger_age, older_age, mid_age)

        # actual initialization
        self.epoch_mid_age_list_young2old \
            = mid_ages_minus_last_epoch \
            + self.epoch_mid_age_list_young2old
        
        # debugging
        # print(self.epoch_mid_age_list_young2old)

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
    def get_feat_by_name(self,
                         feat_name: str,
                         time_idx: ty.Optional[int] = None) \
            -> ty.Tuple[float, ty.List[float]]:

        feat_dict_all_epochs = self.feat_name_epochs_dict[feat_name]

        if time_idx:
            return feat_dict_all_epochs[time_idx]

        else:
            return [v for k, v in feat_dict_all_epochs.items()]
    
    def get_feat_by_type_rel_id(self,
                                feat_type: GeoFeatureType,
                                feat_rel: GeoFeatureRelationship,
                                feat_id: int,
                                time_idx: ty.Optional[int] = None) \
            -> ty.Tuple[float, ty.List[float]]:

        # digging into the 4-nested dictionary
        feat_rel_dict_all_epochs = \
            self.feat_type_rel_featid_epochs_dict[feat_type]
        feat_id_dict_all_epochs = \
            feat_rel_dict_all_epochs[feat_rel]
        feat_dict_all_epochs = \
            feat_id_dict_all_epochs[feat_id]
        
        # print(self.feat_type_rel_featid_epochs_dict.keys())
        # print(feat_rel_dict_all_epochs.keys())
        # print(feat_id_dict_all_epochs.keys())
        # print(feat_dict_all_epochs.keys())
        # print(feat_dict_all_epochs[1])
        
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
    def __call__(self, feat_col: GeoFeatureCollection) -> None:
        ...


class GeoFeatureQuery():
    feat_coll: GeoFeatureCollection
    event_dict: ty.Dict[str,
                        ty.Tuple[ty.List[ty.List[float]],
                                 ty.List[ty.List[ty.List[float]]
                                         ]]]

    def __init__(self, feat_coll) -> None:
        self.feat_coll = feat_coll
        self.event_dict = dict()

    # requirement function 1
    @classmethod
    def cw_feature_equals_value(cls,
                            feat_col: GeoFeatureCollection,
                            val2match: int,
                            feat_name: str = "",
                            feat_id: int = None) \
            -> ty.List[ty.List[float]]:
    
        ft_all_epochs_list: ty.List[GeoFeature] = list()
        
        if feat_name:
            ft_all_epochs_list = feat_col.get_feat_by_name(feat_name)
    
        elif feat_id:
            ft_all_epochs_list = feat_col.get_feat_by_type_rel_id(
                GeoFeatureType.CATEGORICAL,
                GeoFeatureRelationship.WITHIN,
                feat_id
            )

        else:
            exit(("ERROR: Function \'cw_feature_equals_value\' needs a "
                    "categorical-within feature name or feature id."
                    ". Exiting."))

        n_regions = ft_all_epochs_list[0].n_regions  # just grabbing 1st
        mid_ages = feat_col.epoch_mid_age_list_young2old
        
        mid_ages_match_list = list()

        for region_idx in range(n_regions):
            mid_ages_match_this_region = list()
            
            for idx, ft_this_epoch in enumerate(ft_all_epochs_list):
                if ft_this_epoch.value[region_idx] == val2match:
                    mid_ages_match_this_region.append(
                        mid_ages[ft_this_epoch.time_idx-1])

            mid_ages_match_list.append(mid_ages_match_this_region)

        return mid_ages_match_list
    
    # requirement function 2
    @classmethod
    def qw_feature_threshold(cls,
                            feat_col: GeoFeatureCollection,
                            threshold: float,
                            above: bool,
                            feat_name: str = "",
                            feat_id: int = None) \
            -> ty.List[ty.List[float]]:
    
        ft_all_epochs_list: ty.List[GeoFeature] = list()

        if feat_name:
            ft_all_epochs_list = feat_col.get_feat_by_name(feat_name)

        elif feat_id:
            ft_all_epochs_list = feat_col.get_feat_by_type_rel_id(
                GeoFeatureType.QUANTITATIVE,
                GeoFeatureRelationship.WITHIN,
                feat_id
            )

        else:
            exit(("ERROR: Function \'qw_feature_threshold\' needs a "
                    "quantitative-within feature name or feature id."
                    ". Exiting."))
            
        n_regions = ft_all_epochs_list[0].n_regions  # just grabbing 1st
        mid_ages = feat_col.epoch_mid_age_list_young2old
        
        mid_ages_match_list = list()

        for region_idx in range(n_regions):
            mid_ages_match_this_region = list()
            
            for idx, ft_this_epoch in enumerate(ft_all_epochs_list):
                if above and ft_this_epoch.value[region_idx] > threshold:
                    mid_ages_match_this_region.append(
                        mid_ages[ft_this_epoch.time_idx-1])
                    
                elif not above and ft_this_epoch.value[region_idx] < threshold:
                    mid_ages_match_this_region.append(
                        mid_ages[ft_this_epoch.time_idx-1])

            mid_ages_match_list.append(mid_ages_match_this_region)

        return mid_ages_match_list
        
    # requirement function 3
    @classmethod
    def cb_feature_equals_value(cls,
                                feat_col: GeoFeatureCollection,
                                val2match: int,
                                feat_name: str = "",
                                feat_id: int = None
                                ) \
            -> ty.List[ty.List[float]]:
        
        ft_all_epochs_list: ty.List[GeoFeature] = list()

        if feat_name:
            ft_all_epochs_list = feat_col.get_feat_by_name(feat_name)

        elif feat_id:
            ft_all_epochs_list = feat_col.get_feat_by_type_rel_id(
                GeoFeatureType.CATEGORICAL,
                GeoFeatureRelationship.BETWEEN,
                feat_id
            )

        else:
            exit(("ERROR: Function \'qw_feature_threshold\' needs a "
                    "quantitative-within feature name or feature id."
                    ". Exiting."))
            
        n_regions = ft_all_epochs_list[0].n_regions  # just grabbing 1st
        mid_ages = feat_col.epoch_mid_age_list_young2old

        mid_ages_match_2d_list = list()

        for region_idx1 in range(n_regions):
            mid_ages_match_from_region = list()
            
            for region_idx2 in range(n_regions):
                mid_ages_match_to_region = list()        

                for idx, ft_this_epoch in enumerate(ft_all_epochs_list):
                    if ft_this_epoch \
                        .get_val_from_to(region_idx1, region_idx2) \
                            == val2match:
                        mid_ages_match_to_region.append(
                            mid_ages[ft_this_epoch.time_idx-1])

                mid_ages_match_from_region.append(mid_ages_match_to_region)

            mid_ages_match_2d_list.append(mid_ages_match_from_region)

        return mid_ages_match_2d_list
        
    # requirement function 4
    @classmethod
    def qb_feature_threshold(cls,
                             feat_col: GeoFeatureCollection,
                             threshold: float,
                             above: bool,
                             feat_name: str = "",
                             feat_id: int = None
                             ) \
            -> ty.List[ty.List[float]]:
        
        ft_all_epochs_list: ty.List[GeoFeature] = list()

        if feat_name:
            ft_all_epochs_list = feat_col.get_feat_by_name(feat_name)

        elif feat_id:
            ft_all_epochs_list = feat_col.get_feat_by_type_rel_id(
                GeoFeatureType.QUANTITATIVE,
                GeoFeatureRelationship.BETWEEN,
                feat_id
            )

        else:
            exit(("ERROR: Function \'qw_feature_threshold\' needs a "
                    "quantitative-within feature name or feature id."
                    ". Exiting."))
            
        n_regions = ft_all_epochs_list[0].n_regions  # just grabbing 1st
        mid_ages = feat_col.epoch_mid_age_list_young2old

        mid_ages_match_2d_list = list()

        for region_idx1 in range(n_regions):
            mid_ages_match_from_region = list()
            
            for region_idx2 in range(n_regions):
                mid_ages_match_to_region = list()        

                for idx, ft_this_epoch in enumerate(ft_all_epochs_list):
                    if above and ft_this_epoch \
                        .get_val_from_to(region_idx1, region_idx2) \
                            > threshold:
                        mid_ages_match_to_region.append(
                            mid_ages[ft_this_epoch.time_idx-1])
                        
                    elif not above and ft_this_epoch \
                        .get_val_from_to(region_idx1, region_idx2) \
                            < threshold:
                        mid_ages_match_to_region.append(
                            mid_ages[ft_this_epoch.time_idx-1])

                mid_ages_match_from_region.append(mid_ages_match_to_region)

            mid_ages_match_2d_list.append(mid_ages_match_from_region)

        return mid_ages_match_2d_list
                          
    def feed_event_dict_with(self,
                             geo_event_name: str,
                             requirement_fn: MyCallableType) \
            -> None:
        """Populate event dictionary with new event
        
        Args:
            geo_event_name (str): the name one is giving the event
            requirement_fn (function): this function receives a feature
                collection and a list of integers defining the involved
                regions, and returns a list of times
        """

        self.event_dict[geo_event_name] = requirement_fn


if __name__ == "__main__":
    # sys.argv[1]: feature summary file path
    # sys.argv[2]: age summary file path

    fc = GeoFeatureCollection(sys.argv[1],
                              age_summary_fp=sys.argv[2])

    test_basic_stuff = False
    if test_basic_stuff:
        # (1) can do it like this
        # print(fc.get_feat_by_name("qb_1", time_idx=1))
        # feat_list = fc.get_feat_by_name("qb_1")
        # for feat in feat_list:
        #     print(feat)
        
        # (2) or like this
        print(fc.get_feat_by_type_rel_id(
            GeoFeatureType.QUANTITATIVE,
            GeoFeatureRelationship.BETWEEN,
            1,
            time_idx=1))
        feat_list = fc.get_feat_by_type_rel_id(
            GeoFeatureType.QUANTITATIVE,
            GeoFeatureRelationship.BETWEEN,
            1)
        for feat in feat_list:
            print(feat)
        
        # checking getter
        # print(feat_list[0].get_val_from_to(1, 2))
    
    # testing querying
    test_querying = True
    if test_querying:
        fq = GeoFeatureQuery(fc)
        
        requirement_fn1 = \
            GeoFeatureQuery.cw_feature_equals_value(fc, 1, feat_name="cw_1")
        
        # for all regions, gives all times when it happened
        fq.feed_event_dict_with("ancient_sea", requirement_fn1)
        print(fq.event_dict["ancient_sea"])
        # [[5.0], [5.0], [-inf], [-inf]]

        requirement_fn1_1 = \
            GeoFeatureQuery.cw_feature_equals_value(fc, 1, feat_id=1)
        
        # for all regions, gives all times when it happened
        fq.feed_event_dict_with("ancient_sea", requirement_fn1_1)
        print(fq.event_dict["ancient_sea"])
        # [[-inf], [5.0], [5.0], [-inf]]

        # for all regions, gives all times when it happened
        thresh = 25.0
        requirement_fn2 = \
            GeoFeatureQuery.qw_feature_threshold(fc, thresh, True, feat_name="qw_1")
        
        fq.feed_event_dict_with("altitude", requirement_fn2)
        print(fq.event_dict["altitude"])

        # for all pairs of regions, gives all times when it happened
        requirement_fn3 = \
            GeoFeatureQuery.cb_feature_equals_value(fc, 1, feat_name="cb_1")
        
        fq.feed_event_dict_with("land_bridge", requirement_fn3)
        print(fq.event_dict["altitude"])
        # [[[], [5.0, -inf], [5.0], []], [[5.0, -inf], [], [5.0], []], [[], [5.0], [], [5.0]], [[], [], [5.0], []]]

        # for all pairs of regions, gives all times when it happened
        requirement_fn3_1 = \
            GeoFeatureQuery.cb_feature_equals_value(fc, 1, feat_id=1)
        
        fq.feed_event_dict_with("land_bridge", requirement_fn3_1)
        print(fq.event_dict["land_bridge"])
        # [[[], [5.0, -inf], [5.0], []], [[5.0, -inf], [], [5.0], []], [[], [5.0], [], [5.0]], [[], [], [5.0], []]]

        # for all pairs of regions, gives all times when it happened
        thresh = 15.0
        requirement_fn4 = \
            GeoFeatureQuery.qb_feature_threshold(fc, thresh, True, feat_name="qb_1")
        
        fq.feed_event_dict_with("distance", requirement_fn4)
        print(fq.event_dict["distance"])
        # [[[], [-inf], [5.0], [5.0, -inf]], [[-inf], [], [5.0], []], [[5.0], [5.0], [], [5.0]], [[], [], [5.0], []]]
         
        requirement_fn4_1 = \
            GeoFeatureQuery.qb_feature_threshold(fc, thresh, True, feat_id=2)
        
        fq.feed_event_dict_with("distance", requirement_fn4_1)
        print(fq.event_dict["distance"])
        # [[[], [5.0], [-inf], [5.0, -inf]], [[5.0], [], [-inf], []], [[-inf], [-inf], [], [-inf]], [[], [], [-inf], []]]