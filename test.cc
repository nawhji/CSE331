#include <windows.h>
#include <psapi.h>
#include <iostream>
#include <cmath>
#include <ctime>
#include <cstdlib>
using namespace std;

// =============================
void bubble_sort(int A[], int size) {
    PROCESS_MEMORY_COUNTERS pmc_before, pmc_peak;
        
    SIZE_T mem_before = 0, mem_peak = 0;
    cout.flush();
    Sleep(10);
    cout.flush();
    if (GetProcessMemoryInfo(GetCurrentProcess(), &pmc_before, sizeof(pmc_before)))
        mem_before = pmc_before.WorkingSetSize;

    for (int i = 1; i <= size - 1; i++) {
        for (int j = size; j >= i + 1; --j) {
            if (A[j] < A[j - 1]) {
                int temp = A[j - 1];
                A[j - 1] = A[j];
                A[j] = temp;
            }
        }
    }

    
    cout.flush();
    Sleep(10);
    cout.flush();
    if (GetProcessMemoryInfo(GetCurrentProcess(), &pmc_peak, sizeof(pmc_peak)))
        mem_peak = pmc_peak.PeakWorkingSetSize;

    SIZE_T used = (mem_peak > mem_before) ? (mem_peak - mem_before) : 0;
    cout << "[bubble_sort size=" << size << "] Peak: " << used << " bytes\n";
}

// =============================
void swap(int A[], int p, int q) {
    int temp = A[p];
    A[p] = A[q];
    A[q] = temp;
}

void cocktail_shaker_sort(int A[], int size) {
    int l = 1;
    int r = size;
    int it = 1;
    bool swapped = true;

    PROCESS_MEMORY_COUNTERS pmc_before, pmc_peak;
        
    SIZE_T mem_before = 0, mem_peak = 0;
    cout.flush();
    Sleep(10);
    cout.flush();
    if (GetProcessMemoryInfo(GetCurrentProcess(), &pmc_before, sizeof(pmc_before)))
        mem_before = pmc_before.WorkingSetSize;

    while (l != r && swapped) {
        swapped = false;
        if (it % 2 == 1) {
            for (int i = l; i < r; i++) {
                if (A[i] > A[i + 1]) {
                    swap(A, i, i + 1);
                    swapped = true;
                }
            }
            --r;
        } else {
            for (int i = r; i > l; --i) {
                if (A[i] < A[i - 1]) {
                    swap(A, i, i - 1);
                    swapped = true;
                }
            }
            l++;
        }
        it++;
    }

    
    cout.flush();
    Sleep(10);
    cout.flush();
    if (GetProcessMemoryInfo(GetCurrentProcess(), &pmc_peak, sizeof(pmc_peak)))
        mem_peak = pmc_peak.PeakWorkingSetSize;

    SIZE_T used = (mem_peak > mem_before) ? (mem_peak - mem_before) : 0;
    cout << "[cocktail_shaker_sort size=" << size << "] Peak: " << used << " bytes\n";
}

// =============================
void heap_sort(int A[], int size) {
    PROCESS_MEMORY_COUNTERS pmc_before, pmc_peak;
    SIZE_T mem_before = 0, mem_peak = 0;

    auto max_heapify = [&](int i, int heap_size, auto&& max_heapify_ref) -> void {
        int l = i * 2;
        int r = i * 2 + 1;
        int largest = i;

        if (l <= heap_size && A[l] > A[largest]) largest = l;
        if (r <= heap_size && A[r] > A[largest]) largest = r;

        if (largest != i) {
            swap(A[i], A[largest]);
            max_heapify_ref(largest, heap_size, max_heapify_ref);
        }
    };
    

    for (int i = size / 2; i >= 1; --i)
        
    max_heapify(i, size, max_heapify);
    cout.flush();
    Sleep(10);
    cout.flush();
    if (GetProcessMemoryInfo(GetCurrentProcess(), &pmc_before, sizeof(pmc_before)))
        mem_before = pmc_before.WorkingSetSize;

    int heap_size = size;
    for (int i = size; i >= 2; --i) {
        swap(A[1], A[i]);
        --heap_size;
        max_heapify(1, heap_size, max_heapify);
    }

    
    cout.flush();
    Sleep(10);
    cout.flush();
    if (GetProcessMemoryInfo(GetCurrentProcess(), &pmc_peak, sizeof(pmc_peak)))
        mem_peak = pmc_peak.PeakWorkingSetSize;
    SIZE_T used = (mem_peak > mem_before) ? (mem_peak - mem_before) : 0;
    cout << "[heap_sort size=" << size << "] Peak: " << used << " bytes\n";
}

// =============================
void intro_sort(int A[], int size) {
    PROCESS_MEMORY_COUNTERS pmc_before, pmc_peak;
    SIZE_T mem_before = 0, mem_peak = 0;

    const int size_threshold = 16;
    int f = 1;
    int b = size + 1;
    int depth_limit = 2 * floor(log2(b - f));

    auto insertion_sort = [&](int l, int r) {
        for (int j = l + 1; j < r; j++) {
            int key = A[j];
            int i = j - 1;
            while (i >= l && A[i] > key) {
                A[i + 1] = A[i];
                --i;
            }
            A[i + 1] = key;
        }
    };

    auto max_heapify = [&](int f, int b, int i, auto&& self) -> void {
        int l = 2 * i - f + 1;
        int r = 2 * i - f + 2;
        int largest = i;

        if (l < b && A[l] > A[largest]) largest = l;
        if (r < b && A[r] > A[largest]) largest = r;

        if (largest != i) {
            swap(A[i], A[largest]);
            self(f, b, largest, self);
        }
    };

    auto build_max_heap = [&](int f, int b) {
        int n = b - f;
        for (int j = n / 2; j >= 1; --j) {
            int i = f + j - 1;
            max_heapify(f, b, i, max_heapify);
        }
    };

    auto heap_sort = [&](int f, int b) {
        build_max_heap(f, b);
        int heap_end = b;
        int n = b - f;
        for (int i = n; i >= 2; --i) {
            swap(A[f], A[heap_end - 1]);
            --heap_end;
            max_heapify(f, heap_end, f, max_heapify);
        }
    };

    auto median_of_3 = [&](int p, int q, int r) -> int {
        if ((A[p] <= A[q] && A[q] <= A[r]) || (A[r] <= A[q] && A[q] <= A[p])) return q;
        if ((A[q] <= A[p] && A[p] <= A[r]) || (A[r] <= A[p] && A[p] <= A[q])) return p;
        return r;
    };

    auto partition = [&](int p, int r, int pivot) -> int {
        swap(A[pivot], A[r - 1]);
        int x = A[r - 1];
        int i = p - 1;
        for (int j = p; j < r - 1; ++j) {
            if (A[j] <= x) {
                ++i;
                swap(A[i], A[j]);
            }
        }
        swap(A[i + 1], A[r - 1]);
        return i + 1;
    };

    auto introsort_loop = [&](int f, int b, int depth_limit, auto&& self) -> void {
        while (b - f > size_threshold) {
            if (depth_limit == 0) {
                heap_sort(f, b);
                return;
            }
            --depth_limit;
            int pivot = median_of_3(f, f + (b - f) / 2, b - 1);
            int p = partition(f, b, pivot);
            self(f, p, depth_limit, self);
            f = p;
        }
        
    };
    cout.flush();
    Sleep(10);
    cout.flush();
    if (GetProcessMemoryInfo(GetCurrentProcess(), &pmc_before, sizeof(pmc_before)))
        mem_before = pmc_before.WorkingSetSize;

    introsort_loop(f, b, depth_limit, introsort_loop);
    insertion_sort(f, b);

    
    cout.flush();
    Sleep(10);
    cout.flush();
    if (GetProcessMemoryInfo(GetCurrentProcess(), &pmc_peak, sizeof(pmc_peak)))
        mem_peak = pmc_peak.PeakWorkingSetSize;

    SIZE_T used = (mem_peak > mem_before) ? (mem_peak - mem_before) : 0;
    cout << "[intro_sort size=" << size << "] Peak: " << used << " bytes\n";
}

#include <optional>
#include <algorithm>

int shift(optional<int> S[], int p, int r, int l_bound, int r_bound) {
    l_bound = max(p, min(r_bound, l_bound));
    r_bound = max(p, min(r, r_bound));

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

void binary_insert(int ins, optional<int> S[], int k) {
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

void rebalance(optional<int> S[], int end_index) {
    int r = end_index;
    int w = r * 2;
    while (r >= 1) {
        S[w] = S[r];
        S[w - 1] = nullopt;
        --r;
        w = w - 2;
    }
}

void library_sort(int A[], int size) {
    PROCESS_MEMORY_COUNTERS pmc_before, pmc_peak;
    SIZE_T mem_before = 0, mem_peak = 0;
    int S_size = size * 2 + 1;
    optional<int>* S = new optional<int>[S_size];

    cout.flush();
    Sleep(10);
    cout.flush();
    if (GetProcessMemoryInfo(GetCurrentProcess(), &pmc_before, sizeof(pmc_before)))
        mem_before = pmc_before.WorkingSetSize;
    S[0] = 0;
    for (int i = 1; i < S_size; i++) S[i] = nullopt;

    for (int i = 1; i <= (int)ceil(log2(size - 1)); i++) {
        int cur_space = (int)pow(2, i - 1);
        int double_space = (int)pow(2, i);
        rebalance(S, cur_space);

        for (int j = (int)pow(2, i - 1); j <= min(size - 1, (int)pow(2, i)); j++) {
            binary_insert(A[j], S, min(double_space, S_size - 1));
        }
    }

    int index = 1;
    for (int i = 1; i < S_size; i++) {
        if (index > size) break;
        if (S[i].has_value()) {
            A[index++] = S[i].value();
        }
        
    }
    cout.flush();
    Sleep(10);
    cout.flush();
    if (GetProcessMemoryInfo(GetCurrentProcess(), &pmc_peak, sizeof(pmc_peak)))
        mem_peak = pmc_peak.PeakWorkingSetSize;

    SIZE_T used = (mem_peak > mem_before) ? (mem_peak - mem_before) : 0;
    cout << "[library_sort size=" << size << "] Peak: " << used << " bytes\n";

    delete[] S;
}

#include <functional>

void merge_sort(int A[], int size) {
    SIZE_T max_peak_used = 0;

    auto merge = [&](int A[], int p, int q, int r) {
        int n1 = q - p;
        int n2 = r - q;

        PROCESS_MEMORY_COUNTERS pmc_before, pmc_peak;
        SIZE_T mem_before = 0, mem_peak = 0;

        
        cout.flush();
        Sleep(10);
        cout.flush();
        if (GetProcessMemoryInfo(GetCurrentProcess(), &pmc_before, sizeof(pmc_before)))
            mem_before = pmc_before.WorkingSetSize;

        int* L = new int[n1 + 1];
        int* R = new int[n2 + 1];

        
        cout.flush();
        Sleep(10);
        cout.flush();
        if (GetProcessMemoryInfo(GetCurrentProcess(), &pmc_peak, sizeof(pmc_peak)))
            mem_peak = pmc_peak.PeakWorkingSetSize;

        SIZE_T used = (mem_peak > mem_before) ? (mem_peak - mem_before) : 0;
        if (used > max_peak_used) max_peak_used = used;

        for (int i = 0; i < n1; ++i)
            L[i] = A[p + i];
        for (int j = 0; j < n2; ++j)
            R[j] = A[q + j];

        L[n1] = INT_MAX;
        R[n2] = INT_MAX;

        int i = 0, j = 0;
        for (int k = p; k < r; ++k) {
            if (L[i] <= R[j]) {
                A[k] = L[i];
                ++i;
            } else {
                A[k] = R[j];
                ++j;
            }
        }

        delete[] L;
        delete[] R;
    };

    std::function<void(int[], int, int)> merge_sort_rec;
    merge_sort_rec = [&](int A[], int p, int r) {
        if (p < r - 1) {
            int q = (p + r) / 2;
            merge_sort_rec(A, p, q);
            merge_sort_rec(A, q, r);
            merge(A, p, q, r);
        }
    };

    merge_sort_rec(A, 1, size + 1);

    // 출력 한 번만
    cout << "[merge_sort size=" << size << "] Peak: " << max_peak_used << " bytes\n";
}


void quick_sort(int A[], int size) {
    SIZE_T max_peak_used = 0;

    auto partition = [](int A[], int p, int r) {
        int x = A[r - 1]; // pivot
        int i = p - 1;
        for (int j = p; j < r - 1; ++j) {
            if (A[j] <= x) {
                ++i;
                swap(A[i], A[j]);
            }
        }
        swap(A[i + 1], A[r - 1]);
        return i + 1;
    };

    function<void(int[], int, int)> quick_sort_rec;
    quick_sort_rec = [&](int A[], int p, int r) {
        if (p < r - 1) {
            PROCESS_MEMORY_COUNTERS pmc_before, pmc_peak;
            SIZE_T mem_before = 0, mem_peak = 0;

            
            cout.flush();
            Sleep(10);
            cout.flush();
            if (GetProcessMemoryInfo(GetCurrentProcess(), &pmc_before, sizeof(pmc_before)))
                mem_before = pmc_before.WorkingSetSize;

            int q = partition(A, p, r);

            
            cout.flush();
            Sleep(10);
            cout.flush();
            if (GetProcessMemoryInfo(GetCurrentProcess(), &pmc_peak, sizeof(pmc_peak)))
                mem_peak = pmc_peak.PeakWorkingSetSize;

            SIZE_T used = (mem_peak > mem_before) ? (mem_peak - mem_before) : 0;
            if (used > max_peak_used) max_peak_used = used;

            quick_sort_rec(A, p, q);
            quick_sort_rec(A, q + 1, r);
        }
    };

    quick_sort_rec(A, 1, size + 1);

    cout << "[quick_sort size=" << size << "] Peak: " << max_peak_used << " bytes\n";
}



void selection_sort(int A[], int size) {
    PROCESS_MEMORY_COUNTERS pmc_before, pmc_peak;
        
    SIZE_T mem_before = 0, mem_peak = 0;
    cout.flush();
    Sleep(10);
    cout.flush();
    if (GetProcessMemoryInfo(GetCurrentProcess(), &pmc_before, sizeof(pmc_before)))
        mem_before = pmc_before.WorkingSetSize;

    for (int j = size; j > 0; --j) {
        int max = A[1];
        int max_index = 1;
        for (int i = 1; i <= j; ++i) {
            if (A[i] > max) {
                max = A[i];
                max_index = i;
            }
        }
        swap(A[j], A[max_index]);
    }

    
    cout.flush();
    Sleep(10);
    cout.flush();
    if (GetProcessMemoryInfo(GetCurrentProcess(), &pmc_peak, sizeof(pmc_peak)))
        mem_peak = pmc_peak.PeakWorkingSetSize;

    SIZE_T used = (mem_peak > mem_before) ? (mem_peak - mem_before) : 0;
    cout << "[selection_sort size=" << size << "] Peak: " << used << " bytes\n";
}


struct Node {
    int data;
    int child_index;
};

void tournament_sort(int A[], int size) {
    PROCESS_MEMORY_COUNTERS pmc_before, pmc_peak;
        
    SIZE_T mem_before = 0, mem_peak = 0;
    cout.flush();
    Sleep(10);
    cout.flush();
    if (GetProcessMemoryInfo(GetCurrentProcess(), &pmc_before, sizeof(pmc_before)))
        mem_before = pmc_before.WorkingSetSize;

    int leaf_count = (int)pow(2, ceil(log2(size)));
    int total_nodes = 2 * leaf_count;
    Node* tree = new Node[total_nodes];

    for (int i = 0; i < total_nodes; ++i) {
        tree[i].data = INT_MAX;
        tree[i].child_index = -1;
    }

    for (int i = 0; i < size; ++i)
        tree[leaf_count + i].data = A[i + 1];

    for (int i = leaf_count - 1; i > 0; --i) {
        int left = 2 * i, right = 2 * i + 1;
        if (tree[left].data < tree[right].data) {
            tree[i].data = tree[left].data;
            tree[i].child_index = left;
        } else {
            tree[i].data = tree[right].data;
            tree[i].child_index = right;
        }
    }

    for (int i = 1; i <= size; ++i) {
        A[i] = tree[1].data;
        int idx = 1;
        while (tree[idx].child_index != -1)
            idx = tree[idx].child_index;
        tree[idx].data = INT_MAX;

        idx /= 2;
        while (idx > 0) {
            int left = 2 * idx, right = 2 * idx + 1;
            if (tree[left].data < tree[right].data) {
                tree[idx].data = tree[left].data;
                tree[idx].child_index = left;
            } else {
                tree[idx].data = tree[right].data;
                tree[idx].child_index = right;
            }
            idx /= 2;
        }
    }

    
    cout.flush();
    Sleep(10);
    cout.flush();
    if (GetProcessMemoryInfo(GetCurrentProcess(), &pmc_peak, sizeof(pmc_peak)))
        mem_peak = pmc_peak.PeakWorkingSetSize;
    SIZE_T used = (mem_peak > mem_before) ? (mem_peak - mem_before) : 0;
    cout << "[tournament_sort size=" << size << "] Peak: " << used << " bytes\n";

    delete[] tree;
}

void insertion_sort(int A[], int size) {
    PROCESS_MEMORY_COUNTERS pmc_before, pmc_peak;
        
    SIZE_T mem_before = 0, mem_peak = 0;
    cout.flush();
    Sleep(10);
    cout.flush();
    if (GetProcessMemoryInfo(GetCurrentProcess(), &pmc_before, sizeof(pmc_before)))
        mem_before = pmc_before.WorkingSetSize;
    for (int j = 2; j <= size; j++) {
        int key = A[j];
        int i = j - 1;
        while (i > 0 && A[i] > key) {
            A[i + 1] = A[i];
            --i;
        }
        A[i + 1] = key;
        
    }
    cout.flush();
    Sleep(10);
    cout.flush();
    if (GetProcessMemoryInfo(GetCurrentProcess(), &pmc_peak, sizeof(pmc_peak)))
        mem_peak = pmc_peak.PeakWorkingSetSize;
    SIZE_T used = (mem_peak > mem_before) ? (mem_peak - mem_before) : 0;
    cout << "[insertion_sort size=" << size << "] Peak: " << used << " bytes\n";
}

void comb_sort(int A[], int size) {
    PROCESS_MEMORY_COUNTERS pmc_before, pmc_peak;
        
    SIZE_T mem_before = 0, mem_peak = 0;
    cout.flush();
    Sleep(10);
    cout.flush();
    if (GetProcessMemoryInfo(GetCurrentProcess(), &pmc_before, sizeof(pmc_before)))
        mem_before = pmc_before.WorkingSetSize;
    int gap = size - 2;
    const double shrink = 1.3;
    bool swapped = true;

    while (gap > 1 || swapped) {
        gap = (int)(gap / shrink);
        if (gap < 1) gap = 1;

        swapped = false;
        for (int i = 1; i + gap <= size; ++i) {
            if (A[i] > A[i + gap]) {
                swap(A, i, i + gap);
                swapped = true;
            }
        }
        
    }
    cout.flush();
    Sleep(10);
    cout.flush();
    if (GetProcessMemoryInfo(GetCurrentProcess(), &pmc_peak, sizeof(pmc_peak)))
        mem_peak = pmc_peak.PeakWorkingSetSize;
    SIZE_T used = (mem_peak > mem_before) ? (mem_peak - mem_before) : 0;
    cout << "[comb_sort size=" << size << "] Peak: " << used << " bytes\n";
}

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

void tim_sort(int A[], int size) {
    PROCESS_MEMORY_COUNTERS pmc_before, pmc_peak;
    SIZE_T mem_before = 0, mem_peak = 0;

    
    cout.flush();
    Sleep(10);
    cout.flush();
    if (GetProcessMemoryInfo(GetCurrentProcess(), &pmc_before, sizeof(pmc_before)))
        mem_before = pmc_before.WorkingSetSize;

    auto galloping_merge = [&](int start1, int end1, int start2, int end2) {
        int len1 = end1 - start1 + 1;
        int len2 = end2 - start2 + 1;

        int i = 0, j = 0;
        int count_left = 0, count_right = 0;

        if (len1 <= len2) {
            std::vector<int> temp(A + start1, A + end1 + 1);
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
                    A[start1 + i + j] = temp[i++];
                    count_left++;
                    count_right = 0;
                } else {
                    A[start1 + i + j] = A[start2 + j++];
                    count_right++;
                    count_left = 0;
                }

                if (count_left >= 7) {
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
    };

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
                galloping_merge(A1, A2, B1, B2);
                run_stack[n - 3] = {A1, B2};
                run_stack.erase(run_stack.begin() + (n - 2));
            } else {
                galloping_merge(B1, B2, C1, C2);
                run_stack[n - 2] = {B1, C2};
                run_stack.pop_back();
            }
        }
    }

    while (run_stack.size() > 1) {
        auto [R2s, R2e] = run_stack.back(); run_stack.pop_back();
        auto [R1s, R1e] = run_stack.back(); run_stack.pop_back();
        galloping_merge(R1s, R1e, R2s, R2e);
        run_stack.emplace_back(R1s, R2e);
    }

    
    cout.flush();
    Sleep(10);
    cout.flush();
    if (GetProcessMemoryInfo(GetCurrentProcess(), &pmc_peak, sizeof(pmc_peak)))
        mem_peak = pmc_peak.PeakWorkingSetSize;

    SIZE_T used = (mem_peak > mem_before) ? (mem_peak - mem_before) : 0;
    cout << "[tim_sort size=" << size << "] Peak: " << used << " bytes\n";
}


// =============================
void run_test(void (*sort_func)(int[], int), const char* name) {
    int sizes[] = {10000};
    for (int s = 0; s < 1; s++) {
        int size = sizes[s];
        int* A = new int[size + 1];
        for (int i = 1; i <= size; ++i)
            A[i] = rand() % 101;
        sort_func(A, size);
        delete[] A;
    }
}

int main() {
    srand((unsigned)time(NULL));
    run_test(bubble_sort, "bubble_sort");
    run_test(comb_sort, "comb_sort");
    run_test(cocktail_shaker_sort, "cocktail_shaker_sort");
    run_test(heap_sort, "heap_sort");
    run_test(insertion_sort, "insertion_sort");
    run_test(intro_sort, "intro_sort");
    run_test(library_sort, "library_sort");
    run_test(merge_sort, "merge_sort");
    run_test(quick_sort, "quick_sort");
    run_test(selection_sort, "selection_sort");
    run_test(tim_sort, "tim_sort");
    run_test(tournament_sort, "tournament_sort");
    return 0;
}