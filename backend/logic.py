from bitarray import bitarray
import time
import psutil
import os

class CMine:
    def __init__(self, input_file, minCP, minRF, maxOR=0.2):
        self.input_file = input_file
        self.minCP = minCP
        self.minRF = minRF
        self.maxOR = maxOR

        self.Df = {}  # Graph ID to bitset of frequency items
        self.L = []   # Final list of patterns with coverage
        self._load_transactions()

    def _load_transactions(self):
        with open(self.input_file, "r") as f:
            lines = f.readlines()
        
        self.numGraphs = len(lines)
        max_item = 0

        transactions = []
        for line in lines:
            items = list(map(int, line.strip().split()))
            transactions.append(items)
            max_item = max(max_item, max(items, default=0))

        self.numFreq = max_item + 1

        for idx, items in enumerate(transactions):
            bitset = bitarray(self.numFreq)
            bitset.setall(0)
            for item in items:
                bitset[item] = 1
            self.Df[idx] = bitset

    def _coverage(self, graph_id):
        return self.Df[graph_id].count() / self.numFreq

    def _pattern_coverage(self, pattern):
        result = bitarray(self.numFreq)
        result.setall(0)
        for g in pattern:
            result |= self.Df[g]
        return result

    def _overlap_ratio(self, pattern):
        last = pattern[-1]
        previous = pattern[:-1]
        prev_cov = self._pattern_coverage(previous)
        last_cov = self.Df[last]

        intersection = prev_cov & last_cov
        union = prev_cov | last_cov
        cs = union.count() / self.numFreq
        ov = intersection.count() / last_cov.count()
        return ov, cs

    def _get_freq1_patterns(self):
        freq = [(i, self._coverage(i)) for i in self.Df if self._coverage(i) >= self.minCP]
        freq_sorted = sorted(freq, key=lambda x: x[1], reverse=True)
        return [[i[0]] for i in freq_sorted]

    def _join(self, l1, l2):
        new_patterns = []
        for i in range(len(l1)):
            for j in range(i + 1, len(l2)):
                if l1[i][:-1] == l2[j][:-1]:
                    if self._coverage(l1[i][-1]) >= self._coverage(l2[j][-1]):
                        new_pattern = l1[i] + [l2[j][-1]]
                    else:
                        new_pattern = l2[j] + [l1[i][-1]]

                    ov, cs = self._overlap_ratio(new_pattern)
                    if ov <= self.maxOR:
                        if cs >= self.minRF:
                            self.L.append((new_pattern, cs))
                        else:
                            new_patterns.append(new_pattern)
        return new_patterns

    def mine(self):
        self._start_time = time.time()
        freq1 = self._get_freq1_patterns()
        nol = []
        for g in freq1:
            cov = self._coverage(g[0])
            if cov >= self.minRF:
                self.L.append((g, cov))
            else:
                nol.append(g)

        l = 2
        while nol:
            nol = self._join(nol, nol)
            l += 1

        process = psutil.Process(os.getpid())
        self._memoryUSS = process.memory_full_info().uss
        self._memoryRSS = process.memory_info().rss
        self._end_time = time.time()

        self.L.sort(key=lambda x: (x[1], len(x[0])), reverse=True)

    def get_patterns(self):
        return self.L

    def get_memory_usage(self):
        return {"USS": self._memoryUSS, "RSS": self._memoryRSS}

    def get_runtime(self):
        return self._end_time - self._start_time

    def write_patterns(self, output_file=None):
        if output_file is None:
            output_file = self.input_file.split(".")[0] + "_results.txt"
        with open(output_file, "w") as f:
            for pattern, cov in self.L:
                f.write(f"{pattern}, Coverage: {cov:.4f}\n")

# Example usage:
# cm = CMine("Toydata.txt", 0.1, 0.3, 0.2)
# cm.mine()
# cm.write_patterns()
# print("Patterns:", cm.get_patterns())
# print("Memory (USS, RSS):", cm.get_memory_usage())
# print("Runtime (s):", cm.get_runtime())
