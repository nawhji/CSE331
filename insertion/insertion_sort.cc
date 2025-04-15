/* 1-indexed algorithm! */
/* 1-indexed algorithm! */
/* 1-indexed algorithm! */
/* 1-indexed algorithm! */
/* 1-indexed algorithm! */
/* 1-indexed algorithm! */

void insertion_sort(int A[], int size) {
    for (int j = 2; j <= size; j++) {
        int key = A[j];
        int i = j - 1;
        while (i > 0 && A[i] > key) {
            A[i + 1] = A[i];
            --i;
        }
        A[i + 1] = key;
    }
}