#include <cmath>
#include <optional>
#include <algorithm>

/* 1-indexed algorithm! */
/* 1-indexed algorithm! */
/* 1-indexed algorithm! */
/* 1-indexed algorithm! */
/* 1-indexed algorithm! */
/* 1-indexed algorithm! */

using namespace std;

int shift(std::optional<int> S[], int p, int r, int l_bound, int r_bound) {
    l_bound = std::max(p, std::min(r_bound, l_bound));
    r_bound = std::max(p, std::min(r, r_bound));

    int x = r_bound;
    while (r > x && S[x].has_value()) ++x;

    if (!S[x].has_value()) {
        while (x > r_bound) {
            S[x] = S[x - 1];
            --x;
        }
        return x;
    } else {
        x = l_bound;
        while (p < x && S[x].has_value()) --x;
        if (!S[x].has_value()) {
            while (x < l_bound) {
                S[x] = S[x + 1];
                ++x;
            }
            return x;
        }
    }
    return -1;
}

void binary_insert(int ins, std::optional<int> S[], int k) {
    int L = 1, R = k, insert_pos = -1;

    while (L <= R) {
        int m = (L + R) / 2;

        if (!S[m].has_value()) {
            int left = m - 1, right = m + 1;
            while (left >= L && !S[left].has_value()) --left;
            while (right <= R && !S[right].has_value()) ++right;

            if (left < L && right > R) {
                insert_pos = m;
                break;
            } else if (left >= L && (right > R || S[left].value() > ins)) {
                R = left;
            } else if (right <= R && (left < L || S[right].value() <= ins)) {
                L = right;
            } else {
                insert_pos = m;
                break;
            }
        } else if (S[m] == ins) {
            while (m <= R && S[m].has_value() && S[m].value() == ins) ++m;
            insert_pos = m;
            break;
        } else if (S[m].value() > ins) {
            R = m - 1;
        } else {
            L = m + 1;
        }
    }

    if (insert_pos == -1)
        insert_pos = L;

    if (S[insert_pos].has_value()) {
        int new_pos = shift(S, 1, k, insert_pos - 1, insert_pos);
        if (new_pos != -1) S[new_pos] = ins;
    } else {
        S[insert_pos] = ins;
    }
}

void rebalance(std::optional<int> S[], int end_index) {
    int r = end_index;
    int w = r * 2;
    while (r >= 1) {
        S[w] = S[r];
        S[w - 1] = std::nullopt;
        --r;
        w = w - 2;
    }
}

void library_sort(int A[], int size) {
    int S_size = size * 2 + 1;
    std::optional<int>* S = new std::optional<int>[S_size];
    S[0] = 0;
    for (int i = 1; i < S_size; i++) S[i] = std::nullopt;

    for (int i = 1; i <= (int)ceil(log2(size - 1)); i++) {
        int cur_space = (int)pow(2, i - 1);
        int double_space = (int)pow(2, i);
        rebalance(S, cur_space);

        for (int j = (int)pow(2, i - 1); j <= std::min(size - 1, (int)pow(2, i)); j++) {
            binary_insert(A[j], S, std::min(double_space, S_size - 1));
        }
    }

    int index = 1;
    for (int i = 1; i < S_size; i++) {
        if (index > size) break;
        if (S[i].has_value()) {
            A[index++] = S[i].value();
        }
    }

    delete[] S;
}
