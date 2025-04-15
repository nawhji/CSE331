#include <windows.h>
#include <psapi.h>
#include <iostream>
#include <fstream>
#include <chrono>
#include <cstdlib>
#include <cstring>
#include <cstdio>
#include <ctime>
#include <vector>
#include <algorithm>
#include <utility>

struct Record {
    int key;
    int order;
    Record(int k = 0, int o = 0) : key(k), order(o) {}
};

void binary_insertion_sort(Record A[], int p, int q, int r) {
    for (int i = q; i <= r; i++) {
        Record key = A[i];
        int L = p;
        int R = i - 1;
        while (L <= R) {
            int m = (L + R) / 2;
            if (A[m].key > key.key)
                R = m - 1;
            else if (A[m].key == key.key) {
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

int scan_for_run(Record A[], int p, int q, int minrun) {
    int cur = p;
    int count = 1;
    bool ascending = false;

    if (p + 1 <= q && A[p].key <= A[p + 1].key) ascending = true;

    if (ascending) {
        while (cur + 1 <= q && A[cur].key <= A[cur + 1].key) {
            count++; cur++;
        }
    } else {
        while (cur + 1 <= q && A[cur].key > A[cur + 1].key) {
            count++; cur++;
        }
        std::reverse(A + p, A + p + count);
    }

    int end_index = std::min(q, p + std::max(count, minrun) - 1);
    binary_insertion_sort(A, p, p + count, end_index);    

    while (end_index < q && A[end_index].key <= A[end_index + 1].key) {
        end_index++;
    }

    return end_index - p + 1;
}

int gallop(const Record* arr, int p, int len, int key, bool is_left) {
    int l = 0, r = len;
    while (l < r) {
        int m = (l + r) / 2;
        if (is_left) {
            if (arr[p + m].key < key) l = m + 1;
            else r = m;
        } else {
            if (arr[p + m].key <= key) l = m + 1;
            else r = m;
        }   
    }
    return l;
}

void galloping_merge(Record A[], int start1, int end1, int start2, int end2) {
    int len1 = end1 - start1 + 1;
    int len2 = end2 - start2 + 1;

    int i = 0, j = 0;
    int count_left = 0, count_right = 0;

    if (len1 <= len2) {
        std::vector<Record> temp(A + start1, A + end1 + 1);
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

            if (temp[i].key <= A[start2 + j].key) {
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
                int advance = gallop(temp.data(), i, len1 - i, A[start2 + j].key, false);
                for (int x = 0; x < advance; ++x)
                    A[start1 + i + j + x] = temp[i + x];
                i += advance;
                count_left = 0;
            } else if (count_right >= 7) {
                int advance = gallop(A, start2 + j, len2 - j, temp[i].key, true);
                for (int x = 0; x < advance; ++x)
                    A[start1 + i + j + x] = A[start2 + j + x];
                j += advance;
                count_right = 0;
            }
        }
    } else {
        std::vector<Record> temp(A + start2, A + end2 + 1);
    
        int i = end1;
        int j = len2 - 1;
        int k = end2;
    
        while (i >= start1 && j >= 0) {
            if (A[i].key > temp[j].key) {
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

void tim_sort(Record A[], int size) {
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

typedef void (*SortFunction)(Record A[], int size);

struct Result {
    char datatype[20];
    int datasize;
    double time;
    SIZE_T memory;
    bool sorted;
    bool stable;
};

Result run_single_experiment(int datasize, const char* datatype, int run_index, SIZE_T prev_memory[4], SortFunction sort_func) {
    int type_index = 0;
    if (strcmp(datatype, "ascending") == 0) type_index = 0;
    else if (strcmp(datatype, "descending") == 0) type_index = 1;
    else if (strcmp(datatype, "random") == 0) type_index = 2;
    else if (strcmp(datatype, "partially sorted") == 0) type_index = 3;

    Result res;
    strcpy(res.datatype, datatype);
    res.datasize = datasize;
    res.time = 0;
    res.memory = 0;
    res.sorted = false;
    res.stable = false;

    std::vector<Record> original;
    Record* arr = new Record[datasize + 1];

    for (int i = 1; i <= datasize; ++i) {
        int key = 0;
        if (strcmp(datatype, "ascending") == 0) key = i % 101;
        else if (strcmp(datatype, "descending") == 0) key = (datasize - i + 1) % 101;
        else if (strcmp(datatype, "random") == 0) key = rand() % 101;
        else key = i % 101;

        arr[i] = Record(key, i);
        original.emplace_back(key, i);
    }

    if (strcmp(datatype, "partially sorted") == 0) {
        int changes = datasize / 20;
        for (int j = 0; j < changes; ++j) {
            int idx = 1 + rand() % datasize;
            int new_key = rand() % 101;
            arr[idx].key = new_key;
            original[idx - 1].key = new_key;
        }
    }

    PROCESS_MEMORY_COUNTERS pmc;
    SIZE_T baseline_mem = 0;
    if (GetProcessMemoryInfo(GetCurrentProcess(), &pmc, sizeof(pmc)))
        baseline_mem = pmc.WorkingSetSize;

    SIZE_T mem_after_alloc = 0;
    if (GetProcessMemoryInfo(GetCurrentProcess(), &pmc, sizeof(pmc)))
        mem_after_alloc = pmc.WorkingSetSize;
    SIZE_T allocated_mem = (mem_after_alloc > baseline_mem) ? (mem_after_alloc - baseline_mem) : 0;
    Sleep(1);
    if (allocated_mem == 0)
        allocated_mem = (datasize + 1) * sizeof(Record);

    auto start = std::chrono::high_resolution_clock::now();
    sort_func(arr, datasize);
    auto end = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double, std::micro> elapsed = end - start;
    res.time = elapsed.count();

    if (run_index == 0) {
        prev_memory[type_index] = allocated_mem;
        res.memory = allocated_mem;
    } else {
        res.memory = prev_memory[type_index];
    }

    bool sorted_flag = true;
    for (int i = 1; i < datasize; ++i) {
        if (arr[i].key > arr[i + 1].key) {
            std::cout << "Unsorted at: " << i << ", key=" << arr[i].key
            << ", next=" << arr[i + 1].key << std::endl;
            sorted_flag = false;
        }
    }

    std::vector<std::vector<int>> expected_order(101);
    for (int i = 1; i <= datasize; ++i)
        expected_order[original[i - 1].key].push_back(original[i - 1].order);

    std::vector<std::vector<int>> after_order(101);
    for (int i = 1; i <= datasize; ++i)
        after_order[arr[i].key].push_back(arr[i].order);

    bool stable_flag = true;
    for (int k = 0; k <= 100; ++k) {
        if (expected_order[k] != after_order[k]) {
            stable_flag = false;
            break;
        }
    }

    res.sorted = sorted_flag;
    res.stable = stable_flag;

    delete[] arr;
    return res;
}

void runTests(SortFunction sort_func, const char* outputFile) {
    std::ofstream ofs(outputFile);
    ofs << "datatype,datasize,time,memory,sorted,stable\n";
    int sizes[4] = {1000, 10000, 100000, 1000000};
    const char* types_arr[4] = {"ascending", "descending", "random", "partially sorted"};

    for (int s = 0; s < 4; s++) {
        int datasize = sizes[s];
        int rep = 10;
        for (int t = 0; t < 4; t++) {
            const char* dtype = types_arr[t];
            SIZE_T memory_rep[4] = {0, 0, 0, 0};
            for (int r = 0; r < rep; r++) {
                Result res = run_single_experiment(datasize, dtype, r, memory_rep, sort_func);
                ofs << res.datatype << "," << res.datasize << "," << res.time << ","
                    << res.memory << "," << (res.sorted ? "true" : "false") << ","
                    << (res.stable ? "true" : "false") << "\n";
                ofs.flush();
                std::cout << "Run " << r << " finished: " << dtype
                          << ", datasize=" << datasize << std::endl;
            }
        }
    }
    ofs.close();
    std::cout << "Test finished, result is in " << outputFile << "✅." << std::endl;
}

int main() {
    srand((unsigned)time(NULL));
    runTests(tim_sort, "results_tim.csv");
    return 0;
}
