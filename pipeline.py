from datetime import datetime
import exceptions
import pandas as pd
import os

class MotifPipeline:
    def __init__(self):
        self.AME_scoring = "max"
        self.FIMO_thresh = "0.0001"
    
    def time_stamp(self):
        now = datetime.now()
        return(now.strftime("%m/%d/%Y, %H:%M:%S"))

    def write_log(self, message):
        with open(self.output_dir + "log.txt", "a+") as log:
            log.write(self.time_stamp() + ":  " + message + "\n")

    def string_name(self):
        return self.output_dir.split("/")[-2]
    
    def set_num_SV_breakpoints(self):
        file_path = self.output_dir + self.string_name() + "_sv.bed"
        if not os.path.isfile(file_path):
            message = ("Error: the following file path <{path}> "
                    "is missing").format(path=file_path)
            raise exceptions.IncorrectPathError(message)
        else:
            with open(file_path) as f:
                for i, l in enumerate(f):
                    pass
                self.num_SV_breakpoints = i
    
    def set_SV_types(self, SV_types_string):
        for SV_type in SV_types_string.split(","):
            if SV_type not in ["inv", "dup", "tra", "del"]:
                message = ("Error: the following SV type <{SV_type}> is incorrect.").format(SV_type=SV_type)
                raise exceptions.WrongArgumentError(message)
        self.SV_types = SV_types_string.split(",")
    
    def set_motif_path(self, motif_path):
        if not os.path.isfile(motif_path):
            message = ("Error: the following file path <{path}> "
                    "is not correct").format(path=motif_path)
            raise exceptions.IncorrectPathError(message)
        else:
            self.motif_path = motif_path
    
    def set_input_dir(self, input_dir):
        if not os.path.isdir(input_dir):
            message = ("Error: the following directory path <{path}> "
                    "is not correct").format(path=input_dir)
            raise exceptions.IncorrectPathError(message)
        else:
            if input_dir[-1] != "/":
                input_dir += "/"
            self.input_dir = input_dir
    
    def set_sample_attr(self, sample_attr_string):
        dict_attr = {}
        for pair in sample_attr_string.split(","):
            if ":" not in pair:
                message = "Error: please use the correct format for attribute listing"
                raise exceptions.WrongArgumentError(message)
            if pair.split(":")[0] in dict_attr.keys():
                dict_attr[pair.split(":")[0]].append(pair.split(":")[1])
            else:
                dict_attr[pair.split(":")[0]] = [pair.split(":")[1]]
        self.sample_attr = dict_attr

    def set_output_dir(self, output_dir):
        try:
            os.makedirs(output_dir, exist_ok=True)
            if output_dir[-1] != "/":
                output_dir += "/"
            self.output_dir = output_dir
        except:
            message = ("Error: the following directory path <{path}> "
                    "is not correct").format(path=output_dir)
            raise exceptions.IncorrectPathError(message)
    
    def set_rand_sv_ratio(self, rand_sv_ratio):
        try:
            self.rand_sv_ratio = int(rand_sv_ratio)
        except:
            message = ("Error: the following random to SV sequences ratio <{AME_scoring}> "
                    "is not correct").format(rand_sv_ratio=rand_sv_ratio)
            raise exceptions.WrongArgumentError(message)

    def set_AME_scoring(self, AME_scoring):
        if AME_scoring in ["max", "avg"]:
            self.AME_scoring = AME_scoring
        else:
            message = ("Error: the following AME scoring method <{AME_scoring}> "
                    "is not correct").format(AME_scoring=AME_scoring)
            raise exceptions.WrongArgumentError(message)
        
    def set_FIMO_thresh(self, FIMO_thresh):
        try:
            float_thresh = float(FIMO_thresh)
            if float_thresh >= 0 and float_thresh <= 1:
                self.FIMO_thresh = FIMO_thresh
            else:
                message = ("Error: the following value for FIMO threshold <{FIMO_thresh}> "
                    "is not correct").format(FIMO_thresh=FIMO_thresh)
                raise exceptions.WrongArgumentError(message)
        except:
            message = ("Error: the following value for FIMO threshold <{FIMO_thresh}> "
                    "is not correct").format(FIMO_thresh=FIMO_thresh)
            raise exceptions.WrongArgumentError(message)
    
    def set_list_bedpe(self, sample_attr_path):
        if not os.path.isfile(sample_attr_path):
            message = ("Error: the following file path <{path}> "
                    "is not correct").format(path=sample_attr_path)
            raise exceptions.IncorrectPathError(message)
        else:
            df_attr = pd.read_csv(sample_attr_path, sep=",")
            list_bedpe = []
            for group in self.sample_attr.keys():
                for attr in self.sample_attr[group]:
                    list_bedpe += df_attr[df_attr[group] == attr][df_attr.columns[0]].tolist()
            self.list_bedpe = list(set(list_bedpe))
            print(self.list_bedpe)
    
    def write_description(self):
        file_prefix = self.output_dir + self.string_name()
        with open(file_prefix + "_results_summary.txt", "a+") as summary:
            summary.write("Summary of results:\n\n")
            summary.write("SV types: {sv_types}\n".format(sv_types=str(self.SV_types)))
            summary.write("Attributes: {attr}\n".format(attr=str(self.sample_attr)))
            summary.write("SV/random ratio: {sv_rand}\n".format(sv_rand=str(self.rand_sv_ratio)))
            summary.write("AME scoring: {AME_scoring}\n".format(AME_scoring=self.AME_scoring))
            summary.write("FIMO threshold: {FIMO_thresh}\n\n".format(FIMO_thresh=self.FIMO_thresh))
