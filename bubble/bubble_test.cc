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

struct Record {
    int key;
    int order;
    Record(int k = 0, int o = 0) : key(k), order(o) {}
    bool operator==(const Record& other) const {
        return key == other.key && order == other.order;
    }
};

void bubble_sort(Record A[], int size) {
    for (int i = 1; i <= size - 1; i++) {
        for (int j = size; j >= i + 1; --j) {
            if (A[j].key < A[j - 1].key) {
                std::swap(A[j], A[j-1]);
            }
        }
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
            sorted_flag = false;
            break;
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
        int rep = (datasize == 1000000) ? 3 : 10;
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
    runTests(bubble_sort, "results_bubble.csv");
    return 0;
}
