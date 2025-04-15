#include <algorithm>
#include <vector>
#include <utility>

/* 1-indexed algorithm! */
/* 1-indexed algorithm! */
/* 1-indexed algorithm! */
/* 1-indexed algorithm! */
/* 1-indexed algorithm! */
/* 1-indexed algorithm! */

using namespace std;

void binary_insertion_sort(int A[], int p, int q, int r) {
    for (int i = q; i <= r; i++) {
        int key = A[i];
        int L = p;
        int R = i - 1;
        while (L <= R) {
            int m = (L + R) / 2;
            if (A[m] > key)
                R = m - 1;
            else if (A[m] == key) {
                L = m + 1;
                break;
            }
            else
                L = m + 1;
        }
        for (int j = i; j > L; --j) A[j] = A[j - 1];
        A[L] = key;
    }
}

int compute_minrun(int n) {
    int r = 0;
    while (n >= 64) {
        if (n % 2 == 1) r = 1;
        n = n / 2; 
    }
    return n + r;
}

int scan_for_run(int A[], int p, int q, int minrun) {
    int cur = p;
    int count = 1;
    bool ascending = false;

    if (p + 1 <= q && A[p] <= A[p + 1]) ascending = true;

    if (ascending) {
        while (cur + 1 <= q && A[cur] <= A[cur + 1]) {
            count++; cur++;
        }
    } else {
        while (cur + 1 <= q && A[cur] > A[cur + 1]) {
            count++; cur++;
        }
        std::reverse(A + p, A + p + count);
    }

    int end_index = std::min(q, p + std::max(count, minrun) - 1);
    binary_insertion_sort(A, p, p + count, end_index);    

    while (end_index < q && A[end_index] <= A[end_index + 1]) {
        end_index++;
    }

    return end_index - p + 1;
}

int gallop(const int* arr, int p, int len, int key, bool is_left) {
    int l = 0, r = len;
    while (l < r) {
        int m = (l + r) / 2;
        if (is_left) {
            if (arr[p + m] < key) l = m + 1;
            else r = m;
        } else {
            if (arr[p + m] <= key) l = m + 1;
            else r = m;
        }   
    }
    return l;
}

void galloping_merge(int A[], int start1, int end1, int start2, int end2) {
    int len1 = end1 - start1 + 1;
    int len2 = end2 - start2 + 1;

    int i = 0, j = 0;
    int count_left = 0, count_right = 0;

    if (len1 <= len2) {
        std::vector<int> temp(A + start1, A + end1 + 1);
        // merge from front
        while (i < len1 || j < len2) {
            if (i == len1) {
                A[start1 + i + j] = A[start2 + j];
                j++;
                continue;
            }
            if (j == len2) {
                A[start1 + i + j] = temp[i];
                i++;
                continue;
            }

            if (temp[i] <= A[start2 + j]) {
                A[start1 + i + j] = temp[i];
                i++;
                count_left++;
                count_right = 0;
            } else {
                A[start1 + i + j] = A[start2 + j];
                j++;
                count_right++;
                count_left = 0;
            }

            if (count_left >= 7/*min*/) {
                int advance = gallop(temp.data(), i, len1 - i, A[start2 + j], false);
                for (int x = 0; x < advance; ++x)
                    A[start1 + i + j + x] = temp[i + x];
                i += advance;
                count_left = 0;
            } else if (count_right >= 7) {
                int advance = gallop(A, start2 + j, len2 - j, temp[i], true);
                for (int x = 0; x < advance; ++x)
                    A[start1 + i + j + x] = A[start2 + j + x];
                j += advance;
                count_right = 0;
            }
        }
    } else {
        std::vector<int> temp(A + start2, A + end2 + 1);
    
        int i = end1;
        int j = len2 - 1;
        int k = end2;
    
        while (i >= start1 && j >= 0) {
            if (A[i] > temp[j]) {
                A[k--] = A[i--];
            } else {
                A[k--] = temp[j--];
            }
        }
    
        while (j >= 0) {
            A[k--] = temp[j--];
        }
    }
    
}

void tim_sort(int A[], int size) {
    int minrun = compute_minrun(size);
    int start = 1;
    std::vector<std::pair<int, int>> run_stack;

    while (start <= size) {
        int run_size = scan_for_run(A, start, size, minrun);
        run_stack.emplace_back(start, start + run_size - 1);
        start += run_size;

        while (true) {
            int n = run_stack.size();
            if (n < 3) break;

            auto [A1, A2] = run_stack[n - 3];
            auto [B1, B2] = run_stack[n - 2];
            auto [C1, C2] = run_stack[n - 1];

            int lenA = A2 - A1 + 1;
            int lenB = B2 - B1 + 1;
            int lenC = C2 - C1 + 1;

            if (lenC > lenA + lenB && lenB > lenA) break;

            if (lenA < lenC) {
                galloping_merge(A, A1, A2, B1, B2);
                run_stack[n - 3] = {A1, B2};
                run_stack.erase(run_stack.begin() + (n - 2));
            } else {
                galloping_merge(A, B1, B2, C1, C2);
                run_stack[n - 2] = {B1, C2};
                run_stack.pop_back();
            }
        }
    }

    while (run_stack.size() > 1) {
        auto [R2s, R2e] = run_stack.back(); run_stack.pop_back();
        auto [R1s, R1e] = run_stack.back(); run_stack.pop_back();
        galloping_merge(A, R1s, R1e, R2s, R2e);
        run_stack.emplace_back(R1s, R2e);
    }
}