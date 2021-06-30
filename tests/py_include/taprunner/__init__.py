# PyUnit-based test runner with TAP reporting.
#
# Vitold Sedyshev <vit1251@gmail.com>
#
# Version 0.1
#
# Bluntly taken and adapted from https://github.com/vit1251/unittest-tap-reporting
# for ViennaRNA Package by Ronny Lorenz <ronny@tbi.univie.ac.at>
#
# License: LGPL
#

import sys
import time
import unittest


class TAPTestResult(unittest.result.TestResult):
    version = 13

    def __init__(self, stream=None, descriptions=None, plan=None):
        super(TAPTestResult, self).__init__()
        #self.stream = sys.stdout
        self.stream = stream
        self._current = 0
        self.descriptions = descriptions
        self.stream.write("TAP version %d\n" % (self.version, ))
        self.stream.write("%d..%d\n" % (1, plan, ))

    def getDescription(self, test):
        doc_first_line = test.shortDescription()
        if doc_first_line:
            return doc_first_line
        else:
            return str(test)

    def addSuccess(self, test):
        #super(TAPTestResult, self).addSuccess(test)
        #print test
        self.stream.write("ok %d - %s\n" % (self._current+1, self.getDescription(test)))
        self.stream.flush()
        self._current += 1

    def addError(self, test, err):
        (exctype, value, tb) = err
        self.stream.write("not ok %d - %s # %s\n" % (self._current+1, self.getDescription(test), value))
        self.stream.flush()
        self._current += 1

    def addFailure(self, test, err):
        (exctype, value, tb) = err
        self.stream.write("not ok %d - %s # %s\n" % (self._current+1, self.getDescription(test), value))
        self.stream.flush()
        self._current += 1

    def addSkip(self, test, reason):
        self.stream.write("not ok %d - %s # SKIP %s\n" % (self._current+1, self.getDescription(test), reason))
        self.stream.flush()
        self._current += 1


class TAPTestRunner(object):
    """A test runner ia abstract class that displays results in user defined form.

    It prints out the names of tests as they are run, errors as they
    occur, and a summary of the results at the end of the test run.
    """
    resultclass = TAPTestResult

    def __init__(self, stream=sys.stdout, descriptions=True, plan=None,
                 failfast=False, buffer=False, resultclass=None):
        self.stream = stream
        self.descriptions = descriptions
        self.plan = plan
        self.failfast = failfast
        self.buffer = buffer
        if resultclass is not None:
            self.resultclass = resultclass

    def _makeResult(self, test):
        if self.plan is None:
            self.plan = test.countTestCases()
        return self.resultclass(self.stream, self.descriptions, self.plan)

    def run(self, test):
        "Run the given test case or test suite."
        result = self._makeResult(test=test)
        #registerResult(result)
        result.failfast = self.failfast
        result.buffer = self.buffer
        startTime = time.time()
        startTestRun = getattr(result, 'startTestRun', None)
        if startTestRun is not None:
            startTestRun()
        try:
            test(result)
        finally:
            stopTestRun = getattr(result, 'stopTestRun', None)
            if stopTestRun is not None:
                stopTestRun()
        stopTime = time.time()
        timeTaken = stopTime - startTime
        return result
