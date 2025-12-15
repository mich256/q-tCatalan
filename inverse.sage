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
        return [colheights_to_binary(part) for part in parts]

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
            result = RationalDyckPath(colheight_to_binary(riffled))
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

def dyck_poset(n):
    """
    Returns the poset of Dyck paths of semilength n by inclusion.
    """
    youngslattice = posets.YoungsLatticePrincipalOrderIdeal(Partition(list(range(n - 1,0,-1)))).dual()
    cover_relations = youngslattice.cover_relations()
    f = {x: x.to_dyck_word(n) for x in youngslattice}
    new_covers = [[f[x], f[y]] for (x,y) in cover_relations]
    newposet = Poset((list(f.values()), new_covers))
    return newposet

def get_roots_under_dp(dw):
    """
    Returns a list of 0,1-vectors. Each vector corresponds to an element in the root poset of type A that corresponds to a box that lies under the Dyck path given by dw. Each vector in the returned list has 0s in all positions except for a contiguous interval of 1s; if this interval stretches from index i-1 to index j-1, then the vector corresponds to the root e_i - e_j.
    """
    rep = to_column_heights(dw)
    rep = [rep[i] - i for i in range(len(rep))]
    n = len(rep)
    boxes = []
    for i in range(n):
        for j in range(1, rep[i]):
            boxes.append((i + 1, j))
    vectors = []
    for (i, j) in boxes:
        interval_start = i - 1
        interval_end = i + j - 2
        vector = [1 if interval_start <= k <= interval_end else 0 for k in range(n-1)]
        vectors.append(vector)
    return vectors

def sum_vectors(v1, v2):
    """
    Sums two 0-1 vectors component-wise, returning a new vector.
    """
    return [v1[i] + v2[i] for i in range(len(v1))]

def find_missing_vectors(existing_vectors, n):
    """
    Given a list of 0,1 vectors where 1s form contiguous intervals and there is at least one 1 in each vector,
    find all possible vectors of the same length with contiguous 1s that are not in the list.
    """
    # Generate all possible contiguous interval vectors
    all_possible = []

    # Generate all contiguous intervals
    for start in range(n):
        for end in range(start, n):
            vector = [0] * n
            for i in range(start, end + 1):
                vector[i] = 1
            all_possible.append(vector)

    # Convert existing vectors to tuples for faster lookup
    existing_set = set(tuple(v) for v in existing_vectors)

    # Find missing vectors
    missing = []
    for vector in all_possible:
        if tuple(vector) not in existing_set:
            missing.append(vector)

    return missing

def check_valid_root(vector):
    """
    Check if a given 0-1 vector corresponds to a valid root (i.e., is a 0-1 vector with contiguous 1s).
    """
    found_one = False
    found_zero_after_one = False
    for bit in vector:
        if bit != 0 and bit != 1:
            return False
        if bit == 1:
            if found_zero_after_one:
                return False
            found_one = True
        else:
            if found_one:
                found_zero_after_one = True
    return found_one

def check_filtered(chain):
    """
    Check if a given chain of Dyck paths under inclusion satisfies the filtered chain conditions.
    """
    roots_chain = []
    missing_roots_chain = []
    for dw in chain:
        roots_chain.append(get_roots_under_dp(dw))
        missing_roots_chain.append(find_missing_vectors(roots_chain[-1], int(len(dw)/2) - 1))
    for i in range(len(chain)):
        for j in range(i, len(chain)):
            if i + 1 + j + 1 <= len(chain):
                for box_i in roots_chain[i]:
                    for box_j in roots_chain[j]:
                        sum_box = sum_vectors(box_i, box_j)
                        if check_valid_root(sum_box) and sum_box not in roots_chain[i + 1 + j + 1 - 1]:
                            return False
            for box_i in missing_roots_chain[i]:
                    for box_j in missing_roots_chain[j]:
                        sum_box = sum_vectors(box_i, box_j)
                        if check_valid_root(sum_box) and sum_box not in missing_roots_chain[min(i + 1 + j + 1 - 1, len(chain) - 1)]:
                            return False
    return True

def filtered_chains(m, n):
    """
    Returns a list of all filtered m-chains of Dyck paths of semilength n.
    Each chain is a tuple (P1, ..., Pm) of Dyck paths such that:
      - the tuple forms an increasing chain in the Dyck path inclusion poset
      - under the bijection of Dyck paths to order ideals in the root poset of type A, the chain corresponds to a filtered chain of order ideals
    """
    from itertools import product
    poset = dyck_poset(n)

    elements = list(poset)

    chain = []
    lam = elements[0]
    filtered_chains = []

    def construct_multichain(lam):
        nonlocal chain, filtered_chains
        if len(chain) == m:
            if check_filtered(chain):
                filtered_chains.append(chain.copy())
            return
        for el in poset.order_filter([lam]):
            chain.append(el)
            construct_multichain(el)
            chain.pop()

    construct_multichain(lam)
    return filtered_chains

def filtered_chains_generator(m, n):
    """
    Generator version that yields filtered m-chains one by one.
    """
    from itertools import product
    poset = dyck_poset(n)

    elements = list(poset)

    chain = []
    lam = elements[0]

    def construct_multichain(lam):
        nonlocal chain
        if len(chain) == m:
            if check_filtered(chain):
                yield chain.copy()
            return
        for el in poset.order_filter([lam]):
            chain.append(el)
            yield from construct_multichain(el)
            chain.pop()

    yield from construct_multichain(lam)

def area_gluing(tup):
    """
    Given an m-tuple of Dyck paths, returns the area vector of the glued m-Dyck path.
    """
    m = len(tup)
    n = len(tup[0]) // 2
    marea = [0] * n
    for dw in tup:
        marea = sum_vectors(marea, DyckWord(dw).to_area_sequence())
    return marea

def get_area_gluing_map_dict(m, n):
    chains = filtered_chains(m, n)
    di = {}
    for chain in chains:
        mdw = area_to_binary(area_gluing(chain), m=len(chain))
        di["".join(map(str, mdw))] = chain
    return di

def get_area_gluing_pairs_generator(m, n):
    """
    Generator that yields (mdw_string, chain) pairs incrementally.
    This allows processing without storing all chains in memory.
    """
    for chain in filtered_chains_generator(m, n):
        mdw = area_to_binary(area_gluing(chain), m=len(chain))
        mdw_string = "".join(map(str, mdw))
        yield mdw_string, chain

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

def area_to_binary(area, m=1):
    """
    Converts an area sequence of an m-Dyck path to binary representation of an m-Dyck path.
    The area sequence gives the area under the path at each step.
    For an m-Dyck path, the slope is m (m horizontal steps for every 1 vertical step).
    """
    n = len(area)
    # extend with a_{n+1} = 0
    a_ext = area + [0]

    word = []
    for i in range(n):
        # add the i-th up-step
        word.append(1)

        # number of down-steps forced after this up-step:
        d_i = m + a_ext[i] - a_ext[i+1]

        word.extend([0] * d_i)

    return word


def colheights_to_binary(heights):
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

def dyck_path_to_lines(dyck_path):
    """
    Convert a Dyck path's pretty print representation to a list of lines.
    """
    from io import StringIO
    import sys

    # Capture the output of pp()
    old_stdout = sys.stdout
    sys.stdout = captured_output = StringIO()
    DyckWord(dyck_path).pp()
    sys.stdout = old_stdout

    # Process the captured output to get lines
    lines = captured_output.getvalue().rstrip().split('\n')
    return lines

def format_dyck_pairs_side_by_side(dt_paths, chain_paths, n):
    """
    Format pairs of Dyck paths side by side.
    Each pair (dt_path, chain_path) is displayed horizontally next to each other.
    """
    result_lines = []

    for i, (dt_path, chain_path) in enumerate(zip(dt_paths, chain_paths)):
        dt_lines = dyck_path_to_lines(dt_path)
        chain_lines = dyck_path_to_lines(chain_path)

        # Ensure both have the same number of lines
        max_lines = max(len(dt_lines), len(chain_lines))
        while len(dt_lines) < max_lines:
            dt_lines.append(' ' * n)
        while len(chain_lines) < max_lines:
            chain_lines.append(' ' * n)

        # Find the width of each path (pad to n characters)
        dt_width = max(len(line) for line in dt_lines) if dt_lines else n
        dt_width = max(dt_width, n)

        # Add pair label
        if i > 0:
            result_lines.append("")  # Empty line between pairs
        result_lines.append(f"Pair {i+1}:")

        # Combine lines side by side
        for dt_line, chain_line in zip(dt_lines, chain_lines):
            # Pad dt_line to consistent width
            padded_dt_line = dt_line.ljust(dt_width)
            combined_line = padded_dt_line + "   " + chain_line  # 3 spaces between
            result_lines.append(combined_line)

    return result_lines

def print_mismatched_chains(n,m):
    import os

    # Create output file
    os.makedirs('mismatch', exist_ok=True)
    filename = f"mismatch/mismatch_{n}_{m}.txt"

    mismatch_count = 0
    chain_count = 0

    print(f"Starting search for mismatched chains (n={n}, m={m})...")
    print(f"Writing results to: {filename}")

    with open(filename, 'w') as f:
        f.write(f"Mismatched chains for n={n}, m={m}\n")
        f.write("=" * 40 + "\n\n")
        f.flush()

        # Process chains incrementally using generator
        for mdw_string, chain in get_area_gluing_pairs_generator(m, n):
            chain_count += 1

            # Print progress every 50 chains
            if chain_count % 50 == 0:
                print(f"Processed {chain_count} chains, found {mismatch_count} mismatches...")

            binary = list(map(int, list(mdw_string)))
            mdp = RationalDyckPath(binary)
            dt = mdp.split()

            # Check for mismatches
            found_mismatch = False
            for i in range(m):
                if DyckWord(dt[i]).to_area_sequence() != DyckWord(chain[i]).to_area_sequence():
                    found_mismatch = True
                    break

            if found_mismatch:
                mismatch_count += 1
                # Write mismatch to file immediately
                f.write(f"Mismatch #{mismatch_count}:\n")
                f.write(f"{mdw_string}\n")

                # Format Dyck paths side by side
                formatted_pairs = format_dyck_pairs_side_by_side(dt, chain, n)
                for line in formatted_pairs:
                    f.write(line + '\n')

                f.write('\n')  # Extra blank line after each mismatch
                f.flush()  # Force write to disk immediately

                # Print progress to console immediately
                print(f"FOUND MISMATCH #{mismatch_count} at chain {chain_count} (mdw: {mdw_string})")

    print(f"\nSearch complete!")
    print(f"Total chains processed: {chain_count}")
    print(f"Total mismatched chains found: {mismatch_count}")
    print(f"Results written to {filename}")
if __name__ == "__main__":
    pass
