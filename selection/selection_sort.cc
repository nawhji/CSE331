/* 1-indexed algorithm! */
/* 1-indexed algorithm! */
/* 1-indexed algorithm! */
/* 1-indexed algorithm! */
/* 1-indexed algorithm! */
/* 1-indexed algorithm! */

void selection_sort(int A[], int size) {
    for (int j = size; j > 0; --j) {
        int max = A[1];
        int max_index = 1;
        for (int i = 1; i <= j; i++) {
            if (A[i] > max) {
                max = A[i];
                max_index = i;
            }
        }
        int temp = A[j];
        A[j] = max;
        A[max_index] = temp;
    }
}