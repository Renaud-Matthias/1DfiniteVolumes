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
        self.time = self._startTime
        self._iter = 0
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
