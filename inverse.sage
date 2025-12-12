import itertools

class RationalDyckPath:
    def __init__(self, l: list):
        self.dyckword = l
        self.vertical = sum(l)
        self.horizontal = len(l) - sum(l)
        self.slope = self.vertical/self.horizontal
        self.n = self.vertical
        self.m = int(self.horizontal / self.vertical)
        self.validate()

    def __str__(self):
        return self.dyckword.__str__()

    def __repr__(self):
        return self.dyckword.__repr__()

    def validate(self):
        v = 0
        h = 0
        for step in self.dyckword:
            if step == 1:
                v += 1
            elif step == 0:
                h += 1
            else:
                raise Exception('Dyck word can only contain 0 and 1')
            if h * self.slope > v:
                raise Exception('Not a valid rational Dyck path')
        return True

    def bounce_path(self, return_v=False):
        bp = []
        v = [0] * (self.m - 1) # this will store the v_i's
        i = 0
        flag = 'v'
        while len(bp) < len(self.dyckword):
            if flag == 'v':
                v.append(0)
                while i < len(bp):
                    if self.dyckword[i] == 1:
                        v[-1] += 1
                        bp.append(1)
                    i += 1
                while self.dyckword[i] == 1:
                    v[-1] += 1
                    bp.append(1)
                    i += 1
                flag = 'h'
            else:
                l = sum(v[-self.m : ])
                bp += [0] * l
                flag = 'v'
        return RationalDyckPath(bp) if not return_v else (RationalDyckPath(bp), v[self.m - 1:])

    def h_from_v(self, v: list):
        """
        Produces a list of h_i's given a list of v_i's.
        """
        v = [0] * (self.m - 1) + v
        return [sum(v[i- self.m + 1:i + 1]) for i in range(self.m - 1, len(v))]

    def split(self):
        """
        Splits the m-Dyck path into an m-tuple of Dyck paths (according to its bounce path).
        """
        bp, v = self.bounce_path(return_v=True)
        h = self.h_from_v(v)
        parts = [[] for _ in range(self.m)]
        columns = to_column_heights(self.dyckword)
        for i, num in enumerate(h):
            new = columns[sum(h[:i]): sum(h[:i]) + num]
            parts[i % self.m].extend(new)
        return [to_binary(part) for part in parts]


class DyckTuple:
    def __init__(self, tup: tuple):
        self.tup = tuple(RationalDyckPath(dp) for dp in tup)
        self.m = len(tup)
        self.n = len(tup[0]) // 2

    def get_bounce_paths(self):
        return [dp.bounce_path() for dp in self.tup]

    def split(self):
        """
        Splits each Dyck path in the m-tuple according to its bounce path.
        Returns a list of m lists, where each list contains sublists corresponding to the segments defined.
        """
        bounce_paths = self.get_bounce_paths()
        bounce_heights = [to_column_heights(bp.dyckword) for bp in bounce_paths]
        dp_heights = [to_column_heights(dp.dyckword) for dp in self.tup]
        split_tuple = []
        for i in range(0, self.m):
            bp = bounce_heights[i]
            curr_height = 0
            li = []
            for j, height in enumerate(bp):
                if height == curr_height:
                    li[-1].append(dp_heights[i][j])
                else:
                    li.append([dp_heights[i][j]])
                curr_height = height
            split_tuple.append(li)
        return split_tuple

    def glue(self, warning_counters=None):
        """
        Glues the m-tuple of Dyck paths into an m-Dyck path (in the bounce path way).
        """
        split_tuple = self.split()
        riffled = riffle_lists(split_tuple)
        prediction = self.predict(riffled)

        # Check if RationalDyckPath(to_binary(riffled)) raises exception iff prediction is False
        try:
            result = RationalDyckPath(to_binary(riffled))
            if not prediction:
                print(f"WARNING: Expected exception but got valid path. Prediction: {prediction}")
                if warning_counters:
                    warning_counters['false_negative'] += 1
            return result
        except Exception as e:
            if prediction:
                print(f"WARNING: Expected valid path but got exception. Prediction: {prediction}, Exception: {e}")
                if warning_counters:
                    warning_counters['false_positive'] += 1
            raise e

    def predict(self, riffled):
        """
        Predicts whether this m-tuple of Dyck paths can be glued into a valid m-Dyck path.
        """
        if sorted(riffled) == riffled:
            return True
        return False

def riffle_lists(lists):
    """
    Riffles a list of lists by taking one element from each list in order.
    Example: [[[1], [3, 3], [4]], [[2, 2], [3], [4]], [[2, 3], [3], [4]]]
    gives [1, 2, 2, 2, 3, 3, 3, 3, 3, 4, 4, 4]"""
    # Transpose the lists, filling missing values with empty lists
    transposed = itertools.zip_longest(*lists, fillvalue=[])
    # Flatten each group and then flatten the entire result
    return list(itertools.chain.from_iterable(itertools.chain.from_iterable(group) for group in transposed))

def primes(dw):
    """
    Splits a Dyck word into its prime components.
    """
    m = int((len(dw) - sum(dw)) / sum(dw))
    primes = []
    h = 0
    v = 0
    for step in dw:
        if step == 1:
            v += 1
        else:
            h += 1
        if v * m == h:
            primes.append(dw[:h + v])
            dw = dw[h + v:]
            h = 0
            v = 0
    return primes

def to_column_heights(dw):
    """
    Converts a Dyck path from binary representation to column heights.

    The height of column i in a Dyck path is the number of boxes beneath its ith horizontal step when the Dyck path is drawn in a rectangular grid.
    """
    heights = []
    h = 0
    for step in dw:
        if step == 1:
            h += 1
        if step == 0:
            heights.append(h)
    return heights

def to_binary(heights):
    """
    Converts a Dyck path from column heights to binary representation.
    """
    validate_heights(heights)
    dw = [1] * heights[0] + [0]
    for i, height in enumerate(heights[1:]):
        i = i + 1
        dw.extend([1] * (height - heights[i - 1]))
        dw.append(0)
    return dw

def validate_heights(heights):
    """
    Check if a list is a valid column heights representation of a Dyck path.
    """
    if not heights == sorted(heights):
        raise Exception('Heights must be non-decreasing')

def generate_tuples(n, m):
    """
    Generate all m-tuples of Dyck paths of semilength n.
    """
    return list(itertools.product(DyckWords(n), repeat=m))

def find(li, i):
    """
    Find all indices of i in li."""
    return [j for j, x in enumerate(li) if x == i]

def cat(n,m=1):
    """
    Compute the (n,m)-Catalan number.
    """
    from sage.all import binomial
    return binomial((m+1)*n, n) // ((m*n)+1)

def check(n,m):
    """
    Check if the prediction condition holds true for all m-tuples of Dyck paths of semilength n.
    """
    tuples = generate_tuples(n, m)
    valid = 0
    warning_counters = {'false_positive': 0, 'false_negative': 0}

    for t in tuples:
        try:
            dt = DyckTuple(t)
            glue = dt.glue(warning_counters)
            # print(glue)
            valid += 1
        except Exception as e:
            # print(f"{t}: {e}")
            continue

    print(f"Valid tuples: {valid} out of {len(tuples)}")
    print(f"Prediction warnings:")
    print(f"  False positives (predicted valid but got exception): {warning_counters['false_positive']}")
    print(f"  False negatives (predicted exception but got valid): {warning_counters['false_negative']}")
    print(f"  Total prediction errors: {warning_counters['false_positive'] + warning_counters['false_negative']}")

    if valid != cat(n, m):
        print(f"Discrepancy found: valid = {valid}, Catalan = {cat(n, m)}")
        raise Exception("Catalan number mismatch")

if __name__ == "__main__":
    pass
