from mixing import Mixing
from mass_splitter import Splitter

class mass_solver():
    def __init__(self, stream1=None, stream2=None, mass_flow1=None, mass_flow2=None, split_ratio=None):            
        self.mixing = Mixing(stream1 = stream1, stream2=stream2, mass_flow1=mass_flow1, mass_flow2=mass_flow2)
        self.splitter = Splitter(stream=self.mixing.stream_out, mass_flow=self.mixing.mass_flow_out, split_ratio=split_ratio)


    def calculate_case1(self):
        strem3_massflow = self.mixing.stream_out

       