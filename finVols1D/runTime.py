"""
Object to handle time
"""


class runTime:

    def __init__(self, timeDict):
        """
        """
        self._startTime = timeDict["startTime"]
        self._endTime = timeDict["endTime"]
        self._dt = timeDict["dt"]
        self._dtSave = timeDict.get("dtSave")
        if self._dtSave==None:
            self._dtSave = self._endTime
            print("dtSave not defined in runTime, "
                  + "only last time step will be saved")
        elif self._dtSave<self._dt:
            self._dtSave = self._dt
        self.time = self._startTime
        self._iter = 0
        # iteration corresponding to end time
        self._lastIter = round(self._endTime / self._dt)
        self.time_1 = None
        self.time_2 = None


    def _updateTime(self, dt):
        """
        """
        self.time_2 = self.time_1
        self.time_1 = self.time
        self.time += dt
        self._iter += 1


    def loop(self):
        """
        return True if end time has not been reached and update time
        """
        if self.time >= self._endTime:
            return False
        else:
            self._updateTime(self._dt)
            return True


    def __repr__(self):
        """return time informations"""
        return f" iteration: {self._iter}, time: {self.time}, dt: {self._dt}"
