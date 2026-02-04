import itertools
from typing import Union

# Examples
# IN_BOX_VIOLATING_CHAIN = [
NON_FILTERED_CHAIN = [[1,1,0,1,0,1,0,0,], [1,1,1,0,0,0,1,0], [1,1,1,0,0,0,1,0], [1,1,1,0,0,1,0,0]]

class RationalDyckPath:
    def __init__(self, l: Union[list, str]):
        self.dyckword = l if type(l) == list else string_to_chain(l)[0]
        self.vertical = sum(self.dyckword)
        self.horizontal = len(self.dyckword) - sum(self.dyckword)
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

    def areaseq(self):
        return binary_to_areaseq(self.dyckword)

def binary_to_areaseq(bin):
    """
    Converts the binary representation of the m-Dyck path to its area sequence.
    """
    if not isinstance(bin, RationalDyckPath):
        bin = RationalDyckPath(bin)
    area_seq = []
    h = 0
    v = 0
    for step in bin.dyckword:
        if step == 1:
            area_seq.append((int(v * bin.m) - h))
            v += 1
        else:
            h += 1
    return area_seq

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
            result = RationalDyckPath(colheights_to_binary(riffled))
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

def string_to_chain(s):
    """
    Converts a string representation of a chain of Dyck paths into a list of Dyck paths in 0-1 format.
    Example input: "11010000, 11100000, 11101000"
    """
    parts = s.split(',')
    chain = []
    for part in parts:
        dw_str = part.strip()
        dw = [int(c) for c in dw_str]
        chain.append(dw)
    return chain

def mdp_to_mdt(mdp):
    """
    Converts an m-Dyck path in 0-1 format into an m-tuple of Dyck paths in 0-1 format.
    """
    return RationalDyckPath(string_to_chain(mdp)[0]).split()

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

def get_boxes_under_dp(dw):
    """
    Returns a list of boxes (i,j) that lie under the Dyck path given by dw.
    Each box is represented as a tuple (i,j) where i is the column index (1-based) and j is the row index (1-based).
    """
    rep = to_column_heights(dw)
    rep = [rep[i] - i - 1 for i in range(len(rep))]
    n = len(rep)
    boxes = []
    for i in range(n):
        for j in range(1, rep[i] + 1):
            boxes.append((i + 1, j + i + 1))
    return boxes

def dp_from_boxes(boxes, n):
    """
    Reconstruct Dyck path given boxes under the Dyck path and return in 0-1 format.
    Note that this function assumes that the boxes correspond to a valid Dyck path.
    """
    col_heights = [i for i in range(1, n + 1)]
    for (i, j) in boxes:
        col_heights[i - 1] = max(col_heights[i - 1], j)
    return colheights_to_binary(col_heights)

def box_from_root(root):
    """
    Given a root vector of length n-1 with contiguous 1s from index i-1 to j-1, return the corresponding box (i,j+1).
    """
    interval_indices = [index for index, value in enumerate(root) if value == 1]
    if not interval_indices:
        raise Exception("Root vector must have at least one '1'")
    i = interval_indices[0] + 1
    j = interval_indices[-1] + 1
    return (i, j + 1)

def m_catalan(n,m):
    """
    Compute the (n,m)-Catalan number.
    """
    from sage.all import binomial
    return binomial((m+1)*n, n) // ((m*n)+1)

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
    Note that this function assumes that inclusion is already satisfied.
    """
    roots_chain = []
    missing_roots_chain = []
    n = int(len(chain[0])/2) - 1  # Calculate once

    # Pre-compute all roots and missing roots
    for dw in chain:
        roots = get_roots_under_dp(dw)
        roots_chain.append(roots)
        missing_roots_chain.append(find_missing_vectors(roots, n))

    # Convert to sets for faster lookup
    roots_sets = [set(tuple(vec) for vec in roots) for roots in roots_chain]
    missing_sets = [set(tuple(vec) for vec in missing) for missing in missing_roots_chain]

    # Check filtered conditions
    chain_len = len(chain)
    for i in range(chain_len):
        for j in range(i, chain_len):
            target_idx = i + 1 + j + 1 - 1
            if target_idx < chain_len:
                target_roots_set = roots_sets[target_idx]
                # Check in-box condition
                for box_i in roots_chain[i]:
                    for box_j in roots_chain[j]:
                        sum_box = tuple(box_i[k] + box_j[k] for k in range(len(box_i)))
                        if check_valid_root(sum_box) and sum_box not in target_roots_set:
                            return False

            # Check out-box condition
            target_missing_idx = min(target_idx, chain_len - 1) if target_idx < chain_len else chain_len - 1
            target_missing_set = missing_sets[target_missing_idx]

            for box_i in missing_roots_chain[i]:
                for box_j in missing_roots_chain[j]:
                    sum_box = tuple(box_i[k] + box_j[k] for k in range(len(box_i)))
                    if check_valid_root(sum_box) and sum_box not in target_missing_set:
                        return False
    return True

def check_filtered_incremental(partial_chain, roots_cache=None, missing_cache=None):
    """
    Incremental version of check_filtered that can reuse computations from shorter chains.
    """
    if roots_cache is None:
        roots_cache = {}
    if missing_cache is None:
        missing_cache = {}

    chain_len = len(partial_chain)
    if chain_len <= 1:
        return True, roots_cache, missing_cache

    # Only compute for new elements
    roots_chain = []
    missing_roots_chain = []
    n = int(len(partial_chain[0])/2) - 1

    for i, dw in enumerate(partial_chain):
        dw_key = tuple(dw)
        if dw_key in roots_cache:
            roots = roots_cache[dw_key]
            missing = missing_cache[dw_key]
        else:
            roots = get_roots_under_dp(dw)
            missing = find_missing_vectors(roots, n)
            roots_cache[dw_key] = roots
            missing_cache[dw_key] = missing

        roots_chain.append(roots)
        missing_roots_chain.append(missing)

    # Convert to sets for faster lookup
    roots_sets = [set(tuple(vec) for vec in roots) for roots in roots_chain]
    missing_sets = [set(tuple(vec) for vec in missing) for missing in missing_roots_chain]

    # Only check conditions involving the newest element
    new_idx = chain_len - 1

    for i in range(chain_len):
        for j in range(i, chain_len):
            # Skip if neither i nor j is the new element
            if i != new_idx and j != new_idx:
                continue

            target_idx = i + 1 + j + 1 - 1
            if target_idx < chain_len:
                target_roots_set = roots_sets[target_idx]

                for box_i in roots_chain[i]:
                    for box_j in roots_chain[j]:
                        sum_box = tuple(box_i[k] + box_j[k] for k in range(len(box_i)))
                        if check_valid_root(sum_box) and sum_box not in target_roots_set:
                            return False, roots_cache, missing_cache

            target_missing_idx = min(target_idx, chain_len - 1) if target_idx < chain_len else chain_len - 1
            target_missing_set = missing_sets[target_missing_idx]

            for box_i in missing_roots_chain[i]:
                for box_j in missing_roots_chain[j]:
                    sum_box = tuple(box_i[k] + box_j[k] for k in range(len(box_i)))
                    if check_valid_root(sum_box) and sum_box not in target_missing_set:
                        return False, roots_cache, missing_cache

    return True, roots_cache, missing_cache

def check_out_box_filtered(chain, return_witness=False):
    """
    Check if a given tuple of Dyck paths satisfies the second filtered chain condition involving the boxes outside the Dyck paths. (This is the almost same as the function `check_filtered` with the first loop removed)
    """
    roots_chain = []
    missing_roots_chain = []
    for dw in chain:
        roots_chain.append(get_roots_under_dp(dw))
        missing_roots_chain.append(find_missing_vectors(roots_chain[-1], int(len(dw)/2) - 1))
    for i in range(len(chain)):
        for j in range(i, len(chain)):
            for box_i in missing_roots_chain[i]:
                    for box_j in missing_roots_chain[j]:
                        sum_box = sum_vectors(box_i, box_j)
                        if check_valid_root(sum_box) and sum_box not in missing_roots_chain[min(i + 1 + j + 1 - 1, len(chain) - 1)]:
                            print(i, j, box_i, box_j, sum_box, (i + 1 + j + 1 - 1) % len(chain))
                            return False, (i, j, box_from_root(box_i), box_from_root(box_j)) if return_witness else False
    return (True, None) if return_witness else True

def check_in_box_filtered(chain, return_witness=False):
    """
    Check if a given chain of Dyck paths under inclusion satisfies the filtered chain condition involving the boxes inside the Dyck paths. (This is the almost the same as the function `check_filtered` with the second loop removed.)
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
                            return False, (i, j, box_from_root(box_i), box_from_root(box_j)) if return_witness else False
    return (True, None) if return_witness else True

def get_all_in_box_filtered_witnesses(chain):
    roots_chain = []
    missing_roots_chain = []
    for dw in chain:
        roots_chain.append(get_roots_under_dp(dw))
        missing_roots_chain.append(find_missing_vectors(roots_chain[-1], int(len(dw)/2) - 1))
    witnesses = []
    for i in range(len(chain)):
        for j in range(i, len(chain)):
            if i + 1 + j + 1 <= len(chain):
                for box_i in roots_chain[i]:
                    for box_j in roots_chain[j]:
                        # if i == j and box_i[0] >= box_j[0]:
                        #     continue
                        sum_box = sum_vectors(box_i, box_j)
                        if check_valid_root(sum_box) and sum_box not in roots_chain[i + 1 + j + 1 - 1]:
                            witnesses.append((i, j, box_from_root(box_i), box_from_root(box_j)))
    return witnesses

def filtered_chains(m, n):
    """
    Returns a list of all filtered m-chains of Dyck paths of semilength n.
    Each chain is a tuple (P1, ..., Pm) of Dyck paths such that:
      - the tuple forms an increasing chain in the Dyck path inclusion poset
      - under the bijection of Dyck paths to order ideals in the root poset of type A, the chain corresponds to a filtered chain of order ideals
    """
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
    import time

    print(f"Starting filtered_chains_generator for m={m}, n={n}")
    start_time = time.time()

    poset_start = time.time()
    poset = dyck_poset(n)
    poset_time = time.time() - poset_start
    print(f"Built Dyck poset in {poset_time:.2f} seconds")

    elements = list(poset)

    chain = []
    lam = elements[0]

    chains_generated = 0
    last_checkpoint_time = time.time()

    batch_size = 100

    def construct_multichain(lam):
        nonlocal chain, chains_generated, last_checkpoint_time
        if len(chain) == m:
            if check_filtered(chain):
                chains_generated += 1
                if chains_generated % batch_size == 0:
                    current_time = time.time()
                    batch_time = current_time - last_checkpoint_time
                    total_time = current_time - start_time
                    print(f"Generated {chains_generated} filtered chains so far... (Last {batch_size} took {batch_time:.2f}s, Total time: {total_time:.2f}s)")
                    last_checkpoint_time = current_time
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
    chains = filtered_chains_generator(m, n)
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

def generate_m_Dyck_paths(n, m):
    """
    Generate all m-Dyck paths of height n.
    """
    from sage.combinat.tamari_lattices import GeneralizedTamariLattice
    return list(GeneralizedTamariLattice(n * m, m))

def random_m_Dyck_path(n, m):
    """
    Generate a random m-Dyck path of height n using uniform random walk construction.

    Args:
        n: height of the path
        m: slope parameter (m horizontal steps for each vertical step)

    Returns:
        A random m-Dyck path as a binary list
    """
    import random

    # We need n up-steps and m*n down-steps
    total_steps = n * (m + 1)
    up_steps = n
    down_steps = m * n

    path = []
    current_height = 0

    # Build path step by step, respecting the constraint that we never go below the line y = x/m
    for step_num in range(total_steps):
        remaining_steps = total_steps - step_num
        remaining_up = up_steps - sum(path)  # Number of 1s added so far
        remaining_down = down_steps - (step_num - sum(path))  # Number of 0s added so far

        # Calculate current position
        steps_so_far = len(path)
        current_up = sum(path)
        current_down = steps_so_far - current_up

        # Check if we can add an up-step without violating constraints later
        can_go_up = (remaining_up > 0) and (current_down * 1 <= (current_up + 1) * m)

        # Check if we can add a down-step without going below the line
        can_go_down = (remaining_down > 0) and ((current_down + 1) * 1 <= current_up * m)

        # If we must take remaining up steps
        if remaining_down == 0:
            path.append(1)
        # If we must take remaining down steps
        elif remaining_up == 0:
            path.append(0)
        # If both are possible, choose randomly
        elif can_go_up and can_go_down:
            path.append(random.choice([0, 1]))
        # If only up is possible
        elif can_go_up:
            path.append(1)
        # If only down is possible
        elif can_go_down:
            path.append(0)
        else:
            # This shouldn't happen in a well-constructed algorithm
            # Fall back to completing with required steps
            if remaining_up > 0:
                path.append(1)
            else:
                path.append(0)

    return path

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

def print_mismatched_chains(n,m, check_gravity_falls=True):
    import os

    # Create output file
    os.makedirs('slide_mismatch', exist_ok=True)
    filename = f"slide_mismatch/mismatch_{n}_{m}.txt"

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
                if check_gravity_falls:
                    try:
                        if slip_n_slide(dt) != chain:
                            print(f"Slip n slide check failed for mismatch #{mismatch_count} at chain {chain} (mdw: {dt})")
                            f.write("Slip n slide check FAILED!\n")
                    except Exception as e:
                        print(f"Slip n slide check raised exception for mismatch #{mismatch_count} at chain {chain} (mdw: {dt}): {e}")
                        f.write("Slip n slide check RAISED EXCEPTION!\n")

                f.write('\n')  # Extra blank line after each mismatch
                f.flush()  # Force write to disk immediately

                # Print progress to console immediately
                print(f"FOUND MISMATCH #{mismatch_count} at chain {chain_count} (mdw: {mdw_string})")

    print(f"\nSearch complete!")
    print(f"Total chains processed: {chain_count}")
    print(f"Total mismatched chains found: {mismatch_count}")
    print(f"Results written to {filename}")

def print_matched_chains(n,m, check_gravity_falls=True):
    import os

    # Create output file
    os.makedirs('match', exist_ok=True)
    filename = f"match/match_{n}_{m}.txt"

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
                if DyckWord(dt[i]).to_area_sequence() == DyckWord(chain[i]).to_area_sequence():
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
                if check_gravity_falls:
                    if [set(li) for li in gravity_falls([get_boxes_under_dp(dw) for dw in dt])] != [set(li) for li in [get_boxes_under_dp(dw) for dw in chain]]:
                        print(f"Gravity falls check failed for mismatch #{mismatch_count} at chain {chain} (mdw: {dt})")
                        f.write("Gravity falls check FAILED!\n")

                f.write('\n')  # Extra blank line after each mismatch
                f.flush()  # Force write to disk immediately

                # Print progress to console immediately
                print(f"FOUND MISMATCH #{mismatch_count} at chain {chain_count} (mdw: {mdw_string})")

    print(f"\nSearch complete!")
    print(f"Total chains processed: {chain_count}")
    print(f"Total mismatched chains found: {mismatch_count}")
    print(f"Results written to {filename}")


def gravity_falls(tup):
    """
    Stack the Dyck paths represented in tup and move boxes between them by letting them fall down.
    If a tuple exists in list i but not in list j for j > i,
    move it to the largest such j.

    Args:
        tup: List of lists, where the ith list contains tuples representing boxes under the ith Dyck path and above the diagonal OR list of lists representing a chain of Dyck paths in binary format

    Returns:
        The modified original list after letting the boxes "fall" to the bottom.
    """
    n = len(tup)

    if len(tup[0]) != 0 and type(tup[0][0]) is not tuple:
        tup = [get_boxes_under_dp(dw) for dw in tup]
    # Process each list from first to second-to-last
    for i in range(n - 1):
        tuples_to_move = []

        # Check each tuple in the current list
        for tuple_item in tup[i][:]:  # Use slice to avoid modification during iteration
            # Find the largest j > i where this tuple doesn't exist
            target_j = None

            # Check from the end backwards to find the largest valid j
            for j in range(n - 1, i, -1):
                if tuple_item not in tup[j]:
                    target_j = j
                    break

            # If we found a target list, move the tuple there
            if target_j is not None:
                tuples_to_move.append((tuple_item, target_j))

        # Actually move the tuples
        for tuple_item, target_j in tuples_to_move:
            tup[i].remove(tuple_item)
            tup[target_j].append(tuple_item)

    return tup

def slide(chain):
    """
    Given a chain of Dyck paths, perform the gravity falls and then the slide operation on it.
    """
    boxes_list = [get_boxes_under_dp(dw) for dw in chain]
    new_boxes_list = gravity_falls(boxes_list) # it looks like gravity_falls modifies in place, so we might not need to assign it...
    new_chain = [dp_from_boxes(boxes, len(chain[0]) // 2) for boxes in new_boxes_list]
    success, witness = check_in_box_filtered(new_chain, return_witness=True)
    if success:
        return new_chain
    # recall: witness = (i, j, box_i, box_j) where box_i and box_j are (col, row) tuples
    right_box_index = 1 if witness[2][0] <= witness[3][0] else 0
    right_box = witness[2 + right_box_index]
    right_box_dp = witness[right_box_index]
    sliding_boxes_list = [boxes.copy() for boxes in new_boxes_list]


    while gravity_falls(sliding_boxes_list) == new_boxes_list:
        stack_to_slide = [box for box in new_boxes_list[right_box_dp] if box[0] == right_box[0] and box[1]  >= right_box[1]]
        for box in stack_to_slide:
            sliding_boxes_list[right_box_dp].remove(box)
            i = 1
            while (box[0] - i, box[1]) in sliding_boxes_list[right_box_dp + 1]:
                i += 1
            sliding_boxes_list[right_box_dp + 1].append((box[0] - i, box[1]))  # slide left and drop down one level
    new_chain = [dp_from_boxes(boxes, len(chain[0]) // 2) for boxes in sliding_boxes_list]

    return new_chain

def slide_hole(chain):
    """
    Given a chain of Dyck paths, slide a hole that violates the out-box filtered chain condition.
    """
    boxes_list = [get_boxes_under_dp(dw) for dw in chain]
    new_boxes_list = gravity_falls(boxes_list) # this might not be needed, but we keep it for safety
    missing_boxes_list = [[box_from_root(vec) for vec in find_missing_vectors(get_roots_under_dp(dp))] for dp in chain]
    new_chain = [dp_from_boxes(boxes, len(chain[0]) // 2) for boxes in new_boxes_list]
    success, witness = check_out_box_filtered(new_chain, return_witness=True)
    if success:
        return new_chain
    # recall: witness = (i, j, box_i, box_j) where box_i and box_j are (col, row) tuples
    right_box_index = 1 if witness[2][0] <= witness[3][0] else 0
    right_box = witness[2 + right_box_index]
    right_box_dp = witness[right_box_index]
    sliding_boxes_list = [boxes.copy() for boxes in new_boxes_list]

    while gravity_falls(sliding_boxes_list) == new_boxes_list:
        stack_to_slide = [box for box in new_boxes_list[right_box_dp] if box[0] == right_box[0] and box[1] >= right_box[1]]
        for box in stack_to_slide:
            sliding_boxes_list[right_box_dp].remove(box)
            i = 1
            while (box[0] - i, box[1]) in sliding_boxes_list[right_box_dp + 1]:
                i += 1
            sliding_boxes_list[right_box_dp + 1].append((box[0] - i, box[1]))  # slide left and drop down one level
    new_chain = [dp_from_boxes(boxes, len(chain[0]) // 2) for boxes in sliding_boxes_list]

    return new_chain

def print_slide(chain):
    """
    Print the result of performing the slide operation on a given chain of Dyck paths.
    Args:
        chain: A chain (list) of Dyck paths in binary format or an m-Dyck path in binary format.
    """
    if type(chain) is str:
        chain = mdp_to_mdt(chain)
    newchain = slide(chain)
    for dp in newchain:
        DyckWord(dp).pp()
    return newchain

def slip_n_slide(chain):
    """
    Given a chain of Dyck paths, perform the gravity falls and slide operation repeatedly until it stabilizes.
    """
    previous_chain = chain
    while True:
        new_chain = slide(previous_chain)
        if new_chain == previous_chain:
            return new_chain
        previous_chain = new_chain

def print_out_box_violating_bounce_chains(n, m):
    import os

    # Create output file
    os.makedirs('out_box_violations', exist_ok=True)
    filename = f"out_box_violations/out_box_violations_{n}_{m}.txt"

    mismatch_count = 0
    chain_count = 0

    print(f"Writing results to: {filename}")

    with open(filename, 'w') as f:
        f.write(f"Chains violating the out-box filtered chain condition for n={n}, m={m}\n")
        f.write("=" * 40 + "\n\n")
        f.flush()

        # Process chains incrementally using generator
        for mdw_string, chain in get_area_gluing_pairs_generator(m, n):
            chain_count += 1

            # Print progress every 50 chains
            if chain_count % 50 == 0:
                print(f"Processed {chain_count} chains, found {mismatch_count} violations...")

            binary = list(map(int, list(mdw_string)))
            mdp = RationalDyckPath(binary)
            dt = mdp.split()

            if not check_out_box_filtered(chain):
                mismatch_count += 1
                # Write mismatch to file immediately
                f.write(f"Violation #{mismatch_count}:\n")
                f.write(f"{mdw_string}\n")

                # Format Dyck paths side by side
                formatted_pairs = format_dyck_pairs_side_by_side(dt, chain, n)
                for line in formatted_pairs:
                    f.write(line + '\n')

                f.write('\n')  # Extra blank line after each mismatch
                f.flush()  # Force write to disk immediately

                # Print progress to console immediately
                print(f"FOUND VIOLATION #{mismatch_count} at chain {chain_count} (mdw: {mdw_string})")

    print(f"\nSearch complete!")
    print(f"Total chains processed: {chain_count}")
    print(f"Total mismatched chains found: {mismatch_count}")
    print(f"Results written to {filename}")

def print_out_box_violating_chains(n, m, chains):
    import os

    # Create output file
    os.makedirs('out_box_violations_general', exist_ok=True)
    filename = f"out_box_violations_general/out_box_violations_general_{n}_{m}.txt"

    mismatch_count = 0
    chain_count = 0

    print(f"Writing results to: {filename}")

    with open(filename, 'w') as f:
        f.write(f"Chains violating the out-box filtered chain condition for n={n}, m={m}\n")
        f.write("=" * 40 + "\n\n")
        f.flush()

        # Process chains incrementally using generator
        for chain in chains:
            chain_count += 1

            # Print progress every 50 chains
            if chain_count % 50 == 0:
                print(f"Processed {chain_count} chains, found {mismatch_count} violations...")


            if not check_out_box_filtered(chain):
                mismatch_count += 1
                # Write mismatch to file immediately
                f.write(f"Violation #{mismatch_count}:\n")

                # Format Dyck paths side by side
                formatted_pairs = format_dyck_pairs_side_by_side(chain, chain, n)
                for line in formatted_pairs:
                    f.write(line + '\n')

                f.write('\n')  # Extra blank line after each mismatch
                f.flush()  # Force write to disk immediately

                # Print progress to console immediately
                print(f"FOUND VIOLATION #{mismatch_count} at chain {chain_count}")

    print(f"\nSearch complete!")
    print(f"Total chains processed: {chain_count}")
    print(f"Total mismatched chains found: {mismatch_count}")
    print(f"Results written to {filename}")

def print_random_bounce_chain(n,m):
    random_n_m = random_m_Dyck_path(Integer(n),Integer(m))
    mdp_n = RationalDyckPath(random_n_m)
    split_n = mdp_n.split()
    for dp in split_n:
        DyckWord(dp).pp()
    print("m-Dyck word:")
    print(random_n_m)
    print(binary_to_areaseq(random_n_m))
    print()
    print("bounce chain:")
    for dp in split_n:
        print(DyckWord(dp).to_area_sequence())
    print()
    for dp in split_n:
        print(dp)

def generate_integer_partitions(n, m):
    """
    Generate all ways to partition integer n into m non-negative increasing parts.
    Returns a list of tuples where each tuple has m elements summing to n.
    """
    def pad(li, target_length):
        return li + [0] * (target_length - len(li))
    partitions = []
    for partition in Partitions(n, max_length=m):
        padded = pad(list(partition), m)
        partitions.append(list(reversed(padded)))
    return partitions


def is_valid_area_sequence(area_seq):
    """
    Check if an area sequence is valid for a Dyck path:
    - Each consecutive difference is at most 1
    - The sequence starts at 0
    """
    if not area_seq:
        return True

    if area_seq[0] != 0:
        return False
    for i in range(len(area_seq) - 1):
        if area_seq[i+1] - area_seq[i] > 1:
            return False

    # print(f"DEBUG Valid area sequence: {area_seq}")
    return True

def area_gluing_inverse(seq, m):
    """
    Convert a sequence (area or binary) representing an m-Dyck path into the unique filtered chain of Dyck paths.

    Args:
        seq: List of integers representing the area sequence or the binary sequence of an m-Dyck path
        m: Number of Dyck paths in the resulting chain

    Returns:
        The unique valid filtered chain as a list of Dyck paths in binary format, or None if not found
    """
    if seq[0] == 1:
        seq = binary_to_areaseq(seq)
    print(f"Starting with area_seq={seq}, m={m}")
    n = len(seq)
    step = 0

    def generate_partitions(pos, current_partitions):
        """
        Recursively generate all valid partitions of the area sequence.

        Args:
            pos: Current position in area_seq
            current_partitions: List of m lists, each representing partial area sequence
        """
        nonlocal step
        step += 1
        if step % 10000 == 0:
            print(f"  Progress: {pos}/{n}, step {step}")

        if pos == n:
            # Convert area sequences to Dyck paths and check if it's a filtered chain
            # print(f"DEBUG At end, checking partitions: {current_partitions}")
            try:
                chain = []
                for i in range(m):
                    if not is_valid_area_sequence(current_partitions[i]):
                        return None  # Skip invalid area sequences

                    # Convert area sequence to binary Dyck path
                    binary_path = area_to_binary(current_partitions[i])
                    chain.append(binary_path)

                # Check if this chain is filtered
                if check_filtered(chain):
                    return chain

            except Exception as e:
                print(f"Exception during conversion: {e}")
                # Skip invalid conversions
                pass
            return None

        # Generate all partitions of area_seq[pos] into m parts
        current_value = seq[pos]
        # print(f"DEBUG At pos {pos}, partitioning value {current_value}")

        for partition in generate_integer_partitions(current_value, m):
            # Check if adding this partition maintains valid area sequences
            valid = True
            new_partitions = [seq.copy() for seq in current_partitions]

            for i in range(m):
                new_partitions[i].append(partition[i])
                valid = is_valid_area_sequence(new_partitions[i])

            if valid:
                # print(f"DEBUG Trying partition {partition}")
                result = generate_partitions(pos + 1, new_partitions)
                if result is not None:
                    return result

        return None

    # Start with empty area sequences for each of the m Dyck paths
    initial_partitions = [[] for _ in range(m)]
    return generate_partitions(0, initial_partitions)

if __name__ == "__main__":
    pass
