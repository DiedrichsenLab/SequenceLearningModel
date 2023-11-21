class Exp:
    def __init__(self):

        # Model parameters
        self.TN = None
        self.numPress = None
        self.window = None
        self.stimulus = None
        self.RT = None

        # these are specified during the simulation
        self.response = []
        self.decisionTime = []
        self.pressTime = []
        self.stimTime = []
        self.futureTime = []
        
    def reinit(self):
        self.response = []
        self.decisionTime = []
        self.pressTime = []
        self.stimTime = []
        self.futureTime = []

    def singleResp(self):

        self.TN = 1
        self.numPress = 1
        self.window = 1
        self.stimulus = [1]

    def seqShow(self,RT=None):

        self.TN = 1
        self.numPress = 10
        self.window = 10
        self.stimulus = [1,2,3,4,5,1,2,3,4,5]
        if RT is None:
            self.RT = 600
        else:
            self.RT = RT

    def horizon(self, w):

        self.TN = 1
        self.numPress = 13
        self.window = w
        self.stimulus = [1,2,3,4,5,1,2,3,4,5,1,2,3]


