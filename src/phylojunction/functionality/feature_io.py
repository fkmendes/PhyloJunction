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
    str_representation: str

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
        #     self.epoch_age_end_list_young2old
        #     self.epoch_mid_age_list_young2old
        #     self.epoch_age_end_list_old2young
        #     self.epoch_mid_age_list_old2young
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
        self.epoch_age_end_list_old2young \
            = self.epoch_age_end_list_young2old[::-1]

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
        
        self.epoch_mid_age_list_old2young \
            = self.epoch_mid_age_list_young2old[::-1]
        
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
    # bit patterns are in old -> young epochs
    #
    # the value is either a list or a nested list depending on
    # whether the feature was between or within
    geo_cond_bit_dict: ty.Dict[str,
                        ty.Union[ty.List[ty.List[str]],
                                 ty.List[str]]]

    # the value is a 2-D or 3-D list of floats, with all ages
    # of geographic condition change, as indicated by a bit in the
    # bit pattern flipping (if not flips no changes)
    geo_cond_change_times_dict: ty.Dict[str,
                                        ty.Union[ty.List[ty.List[float]],
                                                 ty.List[ty.List[ty.List[float]]]]]

    # this member below helps us tell if geographic conditions were
    # always met, e.g., maybe a barrier never "appears" during the
    # considered epochs, but that's because it existed a long time ago
    # and we will want to know that it did; if it's a "component" barrier
    # of the "full" barrier between complex ranges, we need to know if it's
    # already there so we can say the "full" barrier is observed
    geo_oldest_cond_bit_dict: ty.Dict[str,
                                  ty.Union[ty.List[ty.List[str]],
                                           ty.List[str]]]

    def __init__(self, feat_coll) -> None:
        self.feat_coll = feat_coll
        self.geo_cond_bit_dict = dict()
        self.geo_oldest_cond_bit_dict = dict()
        self.geo_cond_change_times_dict = dict()

    # requirement function 1
    @classmethod
    def cw_feature_equals_value(cls,
                            feat_col: GeoFeatureCollection,
                            val2match: int,
                            feat_name: str = "",
                            feat_id: int = None) \
            -> ty.List[str]:
    
        ft_all_epochs_list: ty.List[GeoFeature] = list()
        
        if feat_name:
            # young -> old
            ft_all_epochs_list = feat_col.get_feat_by_name(feat_name)
    
        elif feat_id:
            # young -> old
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
        # mid_ages = feat_col.epoch_mid_age_list_young2old
        
        # mid_ages_match_list = list()
        bits_match_list = list()

        for region_idx in range(n_regions):
            # mid_ages_match_this_region = list()
            bits_match_this_region_str = str()
            
            # young to old
            # for epoch_idx, ft_this_epoch in enumerate(ft_all_epochs_list):
            for ft_this_epoch in ft_all_epochs_list:
                if ft_this_epoch.value[region_idx] == val2match:
                    # mid_ages_match_this_region.append(
                    #     (epoch_idx, mid_ages[ft_this_epoch.time_idx-1]))
                    bits_match_this_region_str += "1"

                else:
                    bits_match_this_region_str += "0"

            # mid_ages_match_list.append(mid_ages_match_this_region)
            bits_match_list.append(bits_match_this_region_str)

        # return mid_ages_match_list
        return bits_match_list
    
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
            # young -> old
            ft_all_epochs_list = feat_col.get_feat_by_name(feat_name)

        elif feat_id:
            # young -> old
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
        # mid_ages = feat_col.epoch_mid_age_list_young2old
        
        # mid_ages_match_list = list()
        bits_match_list = list()

        for region_idx in range(n_regions):
            # mid_ages_match_this_region = list()
            bits_match_this_region_str = str()

            # for epoch_idx, ft_this_epoch in enumerate(ft_all_epochs_list):
            for ft_this_epoch in ft_all_epochs_list:
                if above:
                    if ft_this_epoch.value[region_idx] > threshold:
                        # mid_ages_match_this_region.append(
                        #     (epoch_idx, mid_ages[ft_this_epoch.time_idx-1]))
                        bits_match_this_region_str += "1"

                    else:
                        bits_match_this_region_str += "0"

                else:
                    if ft_this_epoch.value[region_idx] < threshold:
                        # mid_ages_match_this_region.append(
                        #     (epoch_idx, mid_ages[ft_this_epoch.time_idx-1]))
                        bits_match_this_region_str += "1"

                    else:
                        bits_match_this_region_str += "0"

            # mid_ages_match_list.append(mid_ages_match_this_region)
            bits_match_list.append(bits_match_this_region_str)

        return bits_match_list
        # return mid_ages_match_list
        
    # requirement function 3
    @classmethod
    def cb_feature_equals_value(cls,
                                feat_col: GeoFeatureCollection,
                                val2match: int,
                                feat_name: str = "",
                                feat_id: int = None
                                ) \
            -> ty.List[ty.List[str]]:
            # -> ty.List[ty.List[float]]:
        
        ft_all_epochs_list: ty.List[GeoFeature] = list()

        if feat_name:
            # young -> old
            ft_all_epochs_list = feat_col.get_feat_by_name(feat_name)

        elif feat_id:
            # young -> old
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
        # mid_ages = feat_col.epoch_mid_age_list_young2old

        # mid_ages_match_2d_list = list()
        bits_match_2d_list = list()

        for region_idx1 in range(n_regions):
            # mid_ages_match_from_region = list()
            bits_match_from_region = list()
            
            for region_idx2 in range(n_regions):
                # mid_ages_match_to_region = list()     
                bits_match_to_region_str = str()   

                # for epoch_idx, ft_this_epoch in enumerate(ft_all_epochs_list):
                for ft_this_epoch in ft_all_epochs_list:
                    if ft_this_epoch \
                        .get_val_from_to(region_idx1, region_idx2) \
                            == val2match:
                        # mid_ages_match_to_region.append(
                        #     (epoch_idx, mid_ages[ft_this_epoch.time_idx-1]))
                        bits_match_to_region_str += "1"

                    else:
                        bits_match_to_region_str += "0"

                # mid_ages_match_from_region.append(mid_ages_match_to_region)
                bits_match_from_region.append(bits_match_to_region_str)

            # mid_ages_match_2d_list.append(mid_ages_match_from_region)
            bits_match_2d_list.append(bits_match_from_region)

        # return mid_ages_match_2d_list
        return bits_match_2d_list
        
    # requirement function 4
    @classmethod
    def qb_feature_threshold(cls,
                             feat_col: GeoFeatureCollection,
                             threshold: float,
                             above: bool,
                             feat_name: str = "",
                             feat_id: int = None
                             ) \
            -> ty.List[ty.List[str]]:
            # -> ty.List[ty.List[float]]:
        
        ft_all_epochs_list: ty.List[GeoFeature] = list()

        if feat_name:
            # young -> old
            ft_all_epochs_list = feat_col.get_feat_by_name(feat_name)

        elif feat_id:
            # young -> old
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
        # mid_ages = feat_col.epoch_mid_age_list_young2old

        # mid_ages_match_2d_list = list()
        bits_match_2d_list = list()

        for region_idx1 in range(n_regions):
            # mid_ages_match_from_region = list()
            bits_match_from_region = list()
            
            for region_idx2 in range(n_regions):
                # mid_ages_match_to_region = list()
                bits_match_to_region_str = str()

                # for epoch_idx, ft_this_epoch in enumerate(ft_all_epochs_list):
                for ft_this_epoch in ft_all_epochs_list:
                    if above:
                        if ft_this_epoch .get_val_from_to(region_idx1, region_idx2) \
                                > threshold:
                            # second element is a True/False bit (1 = True, 0 = False)
                            # mid_ages_match_to_region.append(
                            #     (epoch_idx, 1, mid_ages[ft_this_epoch.time_idx-1]))
                            bits_match_to_region_str += "1"

                        else:
                            # mid_ages_match_to_region.append(
                            #     (epoch_idx, 0, None))
                            bits_match_to_region_str += "0"

                    else: 
                        if ft_this_epoch.get_val_from_to(region_idx1, region_idx2) \
                                < threshold:
                            # second element is a True/False bit (1 = True, 0 = False)
                            # mid_ages_match_to_region.append(
                            #     (epoch_idx, 1, mid_ages[ft_this_epoch.time_idx-1]))
                            bits_match_to_region_str += "1"

                        else:
                            # mid_ages_match_to_region.append(
                            #     (epoch_idx, 0, None))
                            bits_match_to_region_str += "0"

                # mid_ages_match_from_region.append(mid_ages_match_to_region)
                # reverse so it's old to young
                bits_match_from_region.append(bits_match_to_region_str[::-1])

            # mid_ages_match_2d_list.append(mid_ages_match_from_region)
            bits_match_2d_list.append(bits_match_from_region)

        # return mid_ages_match_2d_list
        return bits_match_2d_list
                          
    def populate_geo_cond_bit_dicts(
            self,
            geo_cond_name: str,
            requirement_fn: MyCallableType) \
                -> None:
        """Populate geographic condition bit dictionaries

        e.g., for 2 regions, 3 epochs {"altitude": ["010", "010"]}
        the 0's represent a geographic condition that was NOT met,
        while the 1's represent those that were met (depending on
        the requirement function provided by the user)
        
        Args:
            geo_cond_name (str): the name one is giving the geographic
                condition (e.g., \"altitude\", \"distance\", \"ancient
                sea\")
            requirement_fn (function): this function receives a feature
                collection and a list of integers defining the involved
                regions, and returns a nested (1 or 2 dimensions) list
                of strings representing bit patterns
        """

        self.geo_cond_bit_dict[geo_cond_name] = requirement_fn

        self._populate_oldest_geo_cond_bit_dict(geo_cond_name)

        self._populate_geo_cond_change_times_dict(geo_cond_name)

    # internal
    def _populate_geo_cond_change_times_dict(
            self,
            geo_cond_name):
        
        mid_ages_old2young = self.feat_coll.epoch_mid_age_list_old2young
        self.geo_cond_change_times_dict[geo_cond_name] = list()

        # either 1d or 2d list, depending on
        # if within or between feature condition
        geo_cond_bits_list = self.geo_cond_bit_dict[geo_cond_name]

        for region1_str_or_list in geo_cond_bits_list:
            region1_times_list = list()

            # between
            if type(region1_str_or_list) == list:
                # old to young bit
                for region2_str in region1_str_or_list:                    
                    region2_times_list = \
                        [mid_ages_old2young[idx+1] for idx, (i, j) \
                           in enumerate(zip(region2_str,
                                            region2_str[1:])) \
                                                if (i, j) == ('0', '1')]
                    
                    region1_times_list.append(region2_times_list)

                self.geo_cond_change_times_dict[geo_cond_name].append(region1_times_list)

            # within
            else:
                # old to young bit
                region1_times_list.append(
                    [mid_ages_old2young[idx+1] for idx, (i, j) \
                     in enumerate(zip(region1_str_or_list,
                                      region1_str_or_list[1:])) \
                                        if (i, j) == ('0', '1')]
                )

                self.geo_cond_change_times_dict[geo_cond_name].append(region1_times_list)

    def _populate_oldest_geo_cond_bit_dict(self, geo_cond_name):

        self.geo_oldest_cond_bit_dict[geo_cond_name] = list()

        # either 1d or 2d list, depending on
        # if within or between feature condition
        geo_cond_bits_list = self.geo_cond_bit_dict[geo_cond_name]

        for region1_str_or_list in geo_cond_bits_list:
            region1_oldest_bit_list = list()

            # between
            if type(region1_str_or_list) == list:
                for region2_str in region1_str_or_list:                    
                    region2_oldest_bit = region2_str[0]

                    region1_oldest_bit_list.append(region2_oldest_bit)

                self.geo_oldest_cond_bit_dict[geo_cond_name].append(region1_oldest_bit_list)

            # within
            else:
                # old to young bit
                region1_oldest_bit_list = region1_str_or_list[0]

                self.geo_oldest_cond_bit_dict[geo_cond_name].append(region1_oldest_bit_list)

    # getters
    def get_geo_condition_change_times(
            self,
            geo_cond_name):

        if not self.geo_cond_change_times_dict[geo_cond_name]:
            self._populate_geo_cond_change_times_dict(geo_cond_name)

        return self.geo_cond_change_times_dict[geo_cond_name]
    
    def get_geo_oldest_condition_bit(
            self,
            geo_cond_name):

        if not self.geo_oldest_cond_bit_dict[geo_cond_name]:
            self._populate_oldest_geo_cond_bit_dict(geo_cond_name)

        return self.geo_oldest_cond_bit_dict[geo_cond_name]


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
        fq.populate_geo_cond_bit_dicts("ancient_sea", requirement_fn1)
        print("\nancient_sea 1:")
        print(" ".join(fq.geo_cond_bit_dict["ancient_sea"]))
        # ancient_sea 1:
        # 101 101 010 010

        print("\nancient_sea times 1:")
        for k in fq.geo_cond_change_times_dict["ancient_sea"]:
            print(*k)
        # ancient_sea times 1:
        # [5.0]
        # [5.0]
        # [15.0]
        # [15.0]

        print("\nancient_sea oldest bit 1:")
        print(" ".join(fq.geo_oldest_cond_bit_dict["ancient_sea"]))
        # ancient_sea oldest bit 1:
        # 1 1 0 0

        requirement_fn1_1 = \
            GeoFeatureQuery.cw_feature_equals_value(fc, 1, feat_id=1)
        
        # for all regions, gives all times when it happened
        fq.populate_geo_cond_bit_dicts("ancient_sea", requirement_fn1_1)
        print("\nancient_sea 2:")
        print(" ".join(fq.geo_cond_bit_dict["ancient_sea"]))
        # ancient_sea 2:
        # 101 101 010 010

        print("\nancient_sea times 2:")
        for k in fq.geo_cond_change_times_dict["ancient_sea"]:
            print(*k)

        print("\nancient_sea oldest bit 2:")
        print(" ".join(fq.geo_oldest_cond_bit_dict["ancient_sea"]))
        # ancient_sea oldest bit 2:
        # 1 1 0 0

        # for all regions, gives all times when it happened
        thresh = 25.0
        requirement_fn2 = \
            GeoFeatureQuery.qw_feature_threshold(fc, thresh, True, feat_name="qw_1")
        
        fq.populate_geo_cond_bit_dicts("altitude", requirement_fn2)
        print("\naltitude:")
        print(" ".join(fq.geo_cond_bit_dict["altitude"]))
        # altitude:
        # 010 010 101 101

        print("\naltitude times:")
        for k in fq.geo_cond_change_times_dict["altitude"]:
            print(*k)

        print("\naltitude oldest bit:")
        print(" ".join(fq.geo_oldest_cond_bit_dict["altitude"]))
        # altitude oldest bit:
        # 0 0 1 1

        # for all pairs of regions, gives all times when it happened
        requirement_fn3 = \
            GeoFeatureQuery.cb_feature_equals_value(fc, 0, feat_name="cb_1")
        
        fq.populate_geo_cond_bit_dicts("land_bridge", requirement_fn3)
        print("\nland_bridge 1:")
        for k in fq.geo_cond_bit_dict["land_bridge"]:
            print(*k)
        # land_bridge 1:
        # 111 000 010 111
        # 000 111 010 111
        # 111 010 111 010
        # 111 111 010 111

        print("\nland_bridge times 1:")
        for k in fq.geo_cond_change_times_dict["land_bridge"]:
            print(*k)
        # land_bridge times 1:
        # [] [] [15.0] []
        # [] [] [15.0] []
        # [] [15.0] [] [15.0]
        # [] [] [15.0] []

        print("\nland_bridge oldest bit 1:")
        for k in fq.geo_oldest_cond_bit_dict["land_bridge"]:
            print(*k)
        # land_bridge oldest bit 1:
        # 1 0 0 1
        # 0 1 0 1
        # 1 0 1 0
        # 1 1 0 1

        # for all pairs of regions, gives all times when it happened
        requirement_fn3_1 = \
            GeoFeatureQuery.cb_feature_equals_value(fc, 0, feat_id=1)
        
        fq.populate_geo_cond_bit_dicts("land_bridge", requirement_fn3_1)
        print("\nland_bridge 2:")
        for k in fq.geo_cond_bit_dict["land_bridge"]:
            print(*k)
        # land_bridge 2:
        # 111 000 010 111
        # 000 111 010 111
        # 111 010 111 010
        # 111 111 010 111

        print("\nland_bridge times 2:")
        for k in fq.geo_cond_change_times_dict["land_bridge"]:
            print(*k)
        # land_bridge times 2
        # [] [] [15.0] []
        # [] [] [15.0] []
        # [] [15.0] [] [15.0]
        # [] [] [15.0] []

        print("\nland_bridge oldest bit 2:")
        for k in fq.geo_oldest_cond_bit_dict["land_bridge"]:
            print(*k)
        # land_bridge oldest bit 2:
        # 1 0 0 1
        # 0 1 0 1
        # 1 0 1 0
        # 1 1 0 1

        # for all pairs of regions, gives all times when it happened
        thresh = 15.0
        requirement_fn4 = \
            GeoFeatureQuery.qb_feature_threshold(fc, thresh, True, feat_name="qb_1")
        
        fq.populate_geo_cond_bit_dicts("distance", requirement_fn4)
        print("\ndistance 1:")
        for k in fq.geo_cond_bit_dict["distance"]:
            print(*k)
        # distance 1:
        # 000 010 101 111
        # 010 000 101 000
        # 101 101 000 101
        # 000 000 101 000

        print("\ndistance barriers times 1:")
        for k in fq.geo_cond_change_times_dict["distance"]:
            print(*k)
        # distance barriers times 1:
        # [] [15.0] [5.0] []
        # [15.0] [] [5.0] []
        # [5.0] [5.0] [] [5.0]
        # [] [] [5.0] []

        print("\ndistance oldest bit 1:")
        for k in fq.geo_oldest_cond_bit_dict["distance"]:
            print(*k)
        # distance oldest bit 1:
        # 0 0 1 1
        # 0 0 1 0
        # 1 1 0 1
        # 0 0 1 0
         
        requirement_fn4_1 = \
            GeoFeatureQuery.qb_feature_threshold(fc, thresh, True, feat_id=2)
        
        fq.populate_geo_cond_bit_dicts("distance", requirement_fn4_1)
        print("\ndistance 2:")
        for k in fq.geo_cond_bit_dict["distance"]:
            print(*k)
        # distance 2:
        # 000 101 010 110
        # 101 000 010 000
        # 010 010 000 010
        # 000 000 010 000
        
        print("\ndistance barriers times 2:")
        for k in fq.geo_cond_change_times_dict["distance"]:
            print(*k)
        # distance barriers times 2:
        # [] [5.0] [15.0] [15.0]
        # [5.0] [] [15.0] []
        # [15.0] [15.0] [] [15.0]
        # [] [] [15.0] []

        print("\ndistance oldest bit 2:")
        for k in fq.geo_oldest_cond_bit_dict["distance"]:
            print(*k)
        # distance oldest bit 2:
        # 0 1 0 0
        # 1 0 0 0
        # 0 0 0 0
        # 0 0 0 0