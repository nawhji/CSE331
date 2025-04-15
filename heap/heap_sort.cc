/* 1-indexed algorithm! */
/* 1-indexed algorithm! */
/* 1-indexed algorithm! */
/* 1-indexed algorithm! */
/* 1-indexed algorithm! */
/* 1-indexed algorithm! */

void max_heapify(int A[], int i, int heap_size) {
    int l = i * 2;
    int r = i * 2 + 1;
    int largest;

    if (l <= heap_size && A[l] > A[i])
        largest = l;
    else
        largest = i;
    if (r <= heap_size && A[r] > A[largest])
        largest = r;

    if (largest != i) {
        int temp = A[i];
        A[i] = A[largest];
        A[largest] = temp;
        max_heapify(A, largest, heap_size);
    }
}

void build_max_heap(int A[], int size) {
    for (int i = size / 2; i >= 1; --i)
        max_heapify(A, i, size);
}

void heap_sort(int A[], int size) {
    build_max_heap(A, size);
    int heap_size = size;
    for (int i = size; i >= 2; --i) {
        int temp = A[1];
        A[1] = A[i];
        A[i] = temp;
        --heap_size;
        max_heapify(A, 1, heap_size);
    }
}
