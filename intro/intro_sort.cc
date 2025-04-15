#include <cmath>
/* 1-indexed algorithm! */
/* 1-indexed algorithm! */
/* 1-indexed algorithm! */
/* 1-indexed algorithm! */
/* 1-indexed algorithm! */
/* 1-indexed algorithm! */

using namespace std;

const int size_threshold = 16;

////
void insertion_sort(int A[], int l, int r) {
    for (int j = l + 1; j < r; j++) {
        int key = A[j];
        int i = j - 1;
        while (i >= l && A[i] > key) {
            A[i + 1] = A[i];
            --i;
        }
        A[i + 1] = key;
    }
}
////
void max_heapify(int A[], int f, int b, int i) {
    int n = b - f;

    int l = 2 * i - f + 1;
    int r = 2 * i - f + 2;
    int largest = i;

    if (l < b && A[l] > A[largest])
        largest = l;
    if (r < b && A[r] > A[largest])
        largest = r;

    if (largest != i) {
        int temp = A[i];
        A[i] = A[largest];
        A[largest] = temp;
        max_heapify(A, f, b, largest);
    }
}

void build_max_heap(int A[], int f, int b) {
    int n = b - f;
    for (int j = n / 2; j >= 1; --j) {
        int i = f + j - 1;
        max_heapify(A, f, b, i);
    }
}

void heap_sort(int A[], int f, int b) {
    build_max_heap(A, f, b);
    int heap_end = b;
    int n = b - f;
    for (int i = n; i >= 2; --i) {
        int temp = A[f];
        A[f] = A[heap_end - 1];
        A[heap_end - 1] = temp;
        --heap_end;
        max_heapify(A, f, heap_end, f);
    }
}

////
int partition(int A[], int p, int r, int pivot) {
    int x = A[pivot]; // pivot
    int temp = A[pivot];
    A[pivot] = A[r - 1];
    A[r - 1] = temp;

    int i = p - 1;
    for (int j = p; j < r - 1; ++j) {
        if (A[j] <= x) {
            i = i + 1;
            temp = A[j];
            A[j] = A[i];
            A[i] = temp;
        }
    }
    temp = A[i + 1]; // exchange with pivot
    A[i + 1] = A[r - 1];
    A[r - 1] = temp;
    return i + 1;
}

int median_of_3(int A[], int p, int q, int r) {
    if ((A[p] <= A[q] && A[q] <= A[r]) || (A[r] <= A[q] && A[q] <= A[p]))
        return q;
    else if ((A[q] <= A[p] && A[p] <= A[r]) || (A[r] <= A[p] && A[p] <= A[q]))
        return p;
    else
        return r;
}
////

void introsort_loop(int A[], int f, int b, int depth_limit) {
    while (b - f > size_threshold) {
        if (depth_limit == 0) {
            heap_sort(A, f, b);
            return;
        }
        depth_limit = depth_limit - 1;
        int p = partition(A, f, b, median_of_3(A, f, f + (b - f) / 2, b - 1));
        introsort_loop(A, 1, p, depth_limit);
        b = p;
    }
}

void intro_sort1(int A[], int f, int b) {
    introsort_loop(A, f, b, 2 * floor(log2(b - f)));
    insertion_sort(A, f, b);
}

void intro_sort(int A[], int size) {
    intro_sort1(A, 1, size + 1);
}