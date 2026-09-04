from mixing import Mixing
from mass_splitter import Splitter


class Model1():
    def __init__(self, stream1 = None, stream2 = None, mass_flow1 = None, mass_flow2 = None, split_ratio = None):

        self.mixing = Mixing(stream1=stream1, mass_flow1=mass_flow1, stream2=stream2, mass_flow2=mass_flow2)
        self.splitter = Splitter(stream=self.mixing.stream_out, mass_flow=self.mixing.mass_flow_out, split_ratio=split_ratio)

        # Mixer 출구, stream3
        self.stream3 = self.mixing.stream_out
        self.mass_flow3 = self.mixing.mass_flow_out

        # Splitter 출구, stream4
        self.stream4 = self.splitter.stream
        self.mass_flow4 = self.splitter.outlet_mass_flow1

        # Splitter 출구, stream5
        self.stream5 = self.splitter.stream
        self.mass_flow5 = self.splitter.outlet_mass_flow2

