#include "s21_matrix.h"

// проверка матрицы на входные знчения
// 1 - все нормально
int check_matrix(matrix_t *A) {
  int status = 0;
  if (A && A->matrix && A->columns > 0 && A->rows > 0) {
    status = 1;
  }
  return status;
}
// проверка значений на корректность
int nan_or_inf(matrix_t *result) {
  int status = OK;
  for (int i = 0; i < result->rows; i++) {
    for (int j = 0; j < result->columns; j++) {
      if (isnan(result->matrix[i][j]) || isinf(result->matrix[i][j])) {
        status = CALCULATION_ERROR;
        break;
      }
    }
  }
  return status;
}

int s21_create_matrix(int rows, int columns, matrix_t *result) {
  int status = 0;
  if (result && columns > 0 && rows > 0) {
    result->columns = columns;
    result->rows = rows;
    result->matrix = (double **)calloc(rows, sizeof(double *));
    if (result->matrix) {  // проверка на null
      for (int i = 0; i < rows; i++) {
        result->matrix[i] = calloc(columns, sizeof(double));
        if (!result->matrix[i]) {
          status = 1;
          break;
        }
      }
    } else {
      status = 1;
    }

  } else {
    status = 1;
  }
  return status;
}
void s21_remove_matrix(matrix_t *A) {
  if (A && A->matrix) {
    for (int i = 0; i < A->rows; i++) {
      if (A->matrix[i]) {
        free(A->matrix[i]);
      }
    }
    free(A->matrix);
    A->matrix = NULL;
    A->columns = 0;
    A->rows = 0;
  }
}

int s21_eq_matrix(matrix_t *A, matrix_t *B) {
  int status = SUCCESS;
  if (check_matrix(A) && check_matrix(B)) {
    if (A->columns == B->columns && A->rows == B->rows) {
      for (int i = 0; i < A->rows && status; i++) {
        for (int j = 0; j < A->columns; j++) {
          if (fabs(A->matrix[i][j] - B->matrix[i][j]) > EPS) {
            status = FAILURE;
            break;
          }
        }
      }
    } else {
      status = FAILURE;
    }
  }
  return status;
}

int s21_sum_matrix(matrix_t *A, matrix_t *B, matrix_t *result) {
  int status = 0;
  if (check_matrix(A) && check_matrix(B) && result) {
    if (A->columns == B->columns && A->rows == B->rows && !nan_or_inf(A) &&
        !nan_or_inf(B)) {
      s21_create_matrix(A->rows, A->columns, result);
      for (int i = 0; i < A->rows; i++) {
        for (int j = 0; j < A->columns; j++) {
          result->matrix[i][j] = A->matrix[i][j] + B->matrix[i][j];
        }
      }
    } else {
      status = 2;
    }
  } else {
    status = 1;
  }
  return status;
}

int s21_sub_matrix(matrix_t *A, matrix_t *B, matrix_t *result) {
  int status = 0;
  if (check_matrix(A) && check_matrix(B) && result) {
    if (A->columns == B->columns && A->rows == B->rows && !nan_or_inf(A) &&
        !nan_or_inf(B)) {
      s21_create_matrix(A->rows, A->columns, result);
      for (int i = 0; i < A->rows; i++) {
        for (int j = 0; j < A->columns; j++) {
          result->matrix[i][j] = A->matrix[i][j] - B->matrix[i][j];
        }
      }
    } else {
      status = 2;
    }
  } else {
    status = 1;
  }
  return status;
}

int s21_mult_number(matrix_t *A, double number, matrix_t *result) {
  int status = 0;
  if (check_matrix(A) && result) {
    if (!nan_or_inf(A) && !isnan(number) && !isinf(number)) {
      s21_create_matrix(A->rows, A->columns, result);
      for (int i = 0; i < A->rows; i++) {
        for (int j = 0; j < A->columns; j++) {
          result->matrix[i][j] = A->matrix[i][j] * number;
        }
      }
    } else {
      status = 2;
    }
  } else {
    status = 1;
  }
  return status;
}

int s21_mult_matrix(matrix_t *A, matrix_t *B, matrix_t *result) {
  int status = OK;
  if (A == NULL || B == NULL || result == NULL || A->rows < 1 ||
      A->columns < 1 || B->rows < 1 || B->columns < 1)
    status = INCORRECT_MATRIX;
  else if (((A->columns - B->rows) != 0) || nan_or_inf(A) || nan_or_inf(B))
    status = CALCULATION_ERROR;  // A->rows != B->columns возможно это условие
                                 // потребуется
  else {
    s21_create_matrix(A->rows, A->columns, result);
    for (int i = 0; i < A->rows; i++) {
      for (int j = 0; j < B->columns; j++) {
        for (int k = 0; k < B->rows; k++) {
          result->matrix[i][j] += A->matrix[i][k] * B->matrix[k][j];
        }
      }
    }
  }
  return status;
}

int s21_transpose(matrix_t *A, matrix_t *result) {
  int status = 0;
  if (result && check_matrix(A)) {
    if (!nan_or_inf(A)) {
      s21_create_matrix(A->columns, A->rows, result);
      for (int j = 0; j < A->columns; j++) {
        for (int i = 0; i < A->rows; i++) {
          result->matrix[j][i] = A->matrix[i][j];
        }
      }
    } else {
      status = 2;
    }
  } else {
    status = 1;
  }
  return status;
}

void minor_of_matrix(int row, int column, matrix_t *matrix, matrix_t *result) {
  int size = matrix->columns - 1;
  s21_create_matrix(size, size, result);
  int minorRow = 0, minorCol = 0;
  for (int i = 0; i < matrix->rows; i++) {
    if (i != row) {
      for (int j = 0; j < matrix->columns; j++) {
        if (j != column) {
          result->matrix[minorRow][minorCol] = matrix->matrix[i][j];
          minorCol++;
        }
      }
      minorRow++;
      minorCol = 0;
    }
  }
}

int s21_determinant(matrix_t *A, double *result) {
  int status = 0;
  if (check_matrix(A)) {
    if (A->columns == A->rows && !nan_or_inf(A)) {
      if (A->columns == 2) {
        *result = A->matrix[0][0] * A->matrix[1][1] -
                  A->matrix[1][0] * A->matrix[0][1];
      } else if (A->columns == 1) {
        *result = A->matrix[0][0];
      } else {
        *result = 0;
        for (int i = 0; i < A->columns; i++) {
          matrix_t temp = {0};  // n-1 матрицы и ее заполнение
          // s21_create_matrix(A->columns - 1, A->columns - 1, &temp);
          minor_of_matrix(0, i, A, &temp);
          double temp_result;  // временная переменная куда мы записываем
                               // определитель n-1 матрицы
          s21_determinant(&temp, &temp_result);  // status но временный
          s21_remove_matrix(&temp);

          *result += temp_result * A->matrix[0][i] * pow((-1), i);
        }
      }

    } else {
      status = 2;
    }
  } else {
    status = 1;
  }
  return status;
}

int s21_calc_complements(matrix_t *A, matrix_t *result) {
  int status = OK;
  if (A == NULL || result == NULL || A->rows < 1 || A->columns < 1)
    status = INCORRECT_MATRIX;
  else if (A->rows != A->columns)
    status = CALCULATION_ERROR;
  else {
    s21_create_matrix(A->rows, A->columns, result);
    if (A->rows == 1)
      result->matrix[0][0] = A->matrix[0][0];
    else {
      for (int i = 0; i < A->rows; i++) {
        for (int j = 0; j < A->columns; j++) {
          matrix_t minor_matrix;
          minor_of_matrix(i, j, A, &minor_matrix);
          double determinant = 0;
          int sign_of_alg_comp = (((i + j) % 2 == 0) ? 1 : -1);
          s21_determinant(&minor_matrix, &determinant);
          result->matrix[i][j] = sign_of_alg_comp * determinant;
          s21_remove_matrix(&minor_matrix);
        }
      }
    }
  }
  return status;
}

int s21_inverse_matrix(matrix_t *A, matrix_t *result) {
  int status = OK;
  if (A == NULL || result == NULL || A->rows < 1 || A->columns < 1)
    status = INCORRECT_MATRIX;
  else if (A->rows != A->columns)
    status = CALCULATION_ERROR;
  else {
    double determinant = 0.000000;
    s21_determinant(A, &determinant);
    if (fabs(determinant) < EPS)
      status = 2;
    else if (A->rows == 1) {
      s21_create_matrix(1, 1, result);
      result->matrix[0][0] = 1 / determinant;
    } else {
      matrix_t new;
      s21_calc_complements(A, result);
      s21_transpose(result, &new);
      s21_remove_matrix(result);
      s21_mult_number(&new, (1.0 / determinant), result);
      s21_remove_matrix(&new);
    }
  }
  return status;
}