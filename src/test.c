#include <check.h>

#include "s21_matrix.h"

#define SUITE_NAME "Unit"
#define TCASE_NAME "Tests"

START_TEST(test_create_matrix) {
  int rows = 3, columns = 2;
  matrix_t A = {0};
  int res = s21_create_matrix(rows, columns, &A);
  ck_assert_int_eq(res, OK);
  ck_assert_int_eq(A.rows, rows);
  ck_assert_int_eq(A.columns, columns);
  if (!res) s21_remove_matrix(&A);

  int rows2 = 0, columns2 = 1;
  int res2 = s21_create_matrix(rows2, columns2, &A);
  ck_assert_int_eq(res2, INCORRECT_MATRIX);
  if (!res2) s21_remove_matrix(&A);
  ck_assert_int_eq(A.rows, 0);
  ck_assert_int_eq(A.columns, 0);

  int rows3 = 10, columns3 = 10;
  int res3 = s21_create_matrix(rows3, columns3, &A);
  ck_assert_int_eq(res3, OK);
  ck_assert_int_eq(A.rows, rows3);
  ck_assert_int_eq(A.columns, columns3);
  if (!res3) s21_remove_matrix(&A);

  int rows4 = 1, columns4 = 5;
  int res4 = s21_create_matrix(rows4, columns4, &A);
  ck_assert_int_eq(res4, OK);
  ck_assert_int_eq(A.rows, rows4);
  ck_assert_int_eq(A.columns, columns4);
  if (!res4) s21_remove_matrix(&A);

  int rows5 = 10, columns5 = -10;
  int res5 = s21_create_matrix(rows5, columns5, &A);
  ck_assert_int_eq(res5, INCORRECT_MATRIX);
  ck_assert_int_eq(A.rows, 0);
  ck_assert_int_eq(A.columns, 0);
  if (!res5) s21_remove_matrix(&A);
}
END_TEST

START_TEST(test_remove_matrix) {
  int rows = 3, columns = 2;
  matrix_t A = {0};
  s21_create_matrix(rows, columns, &A);
  s21_remove_matrix(&A);
  ck_assert_int_eq(A.rows, 0);
  ck_assert_int_eq(A.columns, 0);

  int rows2 = 5, columns2 = 1;
  s21_create_matrix(rows2, columns2, &A);
  s21_remove_matrix(&A);
  ck_assert_int_eq(A.rows, 0);
  ck_assert_int_eq(A.columns, 0);
  ck_assert_int_eq(A.matrix == NULL, 1);

  int rows3 = 6, columns3 = 6;
  s21_create_matrix(rows3, columns3, &A);
  s21_remove_matrix(&A);
  ck_assert_int_eq(A.rows, 0);
  ck_assert_int_eq(A.columns, 0);
  ck_assert_int_eq(A.matrix == NULL, 1);
}
END_TEST

START_TEST(test_eq_martix) {
  int rows = 2, columns = 3;
  matrix_t A = {0}, B = {0};
  s21_create_matrix(rows, columns, &A);
  s21_create_matrix(rows, columns, &B);
  A.matrix[0][0] = 1.0;
  B.matrix[0][0] = 1.0;
  A.matrix[0][1] = 0.0000007;
  B.matrix[0][1] = 0.0000007;
  A.matrix[0][2] = -16.7865443211;
  B.matrix[0][2] = -16.78654432541;
  A.matrix[1][0] = 7.08;
  B.matrix[1][0] = 7.08;
  A.matrix[1][1] = 3.11111111111;
  B.matrix[1][1] = 3.11111113;
  A.matrix[1][2] = -10.0;
  B.matrix[1][2] = -10.0;
  int res = s21_eq_matrix(&A, &B);
  ck_assert_int_eq(res, 1);
  s21_remove_matrix(&A);
  s21_remove_matrix(&B);

  int rows2 = 2, columns2 = 3;
  int rows3 = 3, columns3 = 2;
  s21_create_matrix(rows2, columns2, &A);
  s21_create_matrix(rows3, columns3, &B);
  int res2 = s21_eq_matrix(&A, &B);
  ck_assert_int_eq(res2, 0);
  s21_remove_matrix(&A);
  s21_remove_matrix(&B);

  int rows4 = 2, columns4 = 2;
  s21_create_matrix(rows4, columns4, &A);
  s21_create_matrix(rows4, columns4, &B);
  A.matrix[0][0] = 3;
  B.matrix[0][0] = 3;
  A.matrix[0][1] = 2;
  B.matrix[0][1] = 2;
  A.matrix[1][0] = 5;
  B.matrix[1][0] = 5;
  A.matrix[1][1] = 4;
  B.matrix[1][1] = 4;
  int res3 = s21_eq_matrix(&A, &B);
  ck_assert_int_eq(res3, 1);
  s21_remove_matrix(&A);
  s21_remove_matrix(&B);

  int rows5 = 10, columns5 = 4;
  int rows6 = 4, columns6 = 10;
  s21_create_matrix(rows5, columns5, &A);
  s21_create_matrix(rows6, columns6, &B);
  int res4 = s21_eq_matrix(&A, &B);
  ck_assert_int_eq(res4, 0);
  s21_remove_matrix(&A);
  s21_remove_matrix(&B);

  int rows7 = 3, columns7 = 3;
  s21_create_matrix(rows7, columns7, &A);
  s21_create_matrix(rows7, columns7, &B);
  A.matrix[0][0] = 1.08765498;
  B.matrix[0][0] = 1.087654986342;
  A.matrix[0][1] = 0.0000007;
  B.matrix[0][1] = 0.0000007;
  A.matrix[0][2] = 3.11111111;
  B.matrix[0][2] = 3.11111111;
  A.matrix[1][0] = 3.8;
  B.matrix[1][0] = 3.8;
  A.matrix[1][1] = -0.01;
  B.matrix[1][1] = -0.01;
  A.matrix[1][2] = 875421.09875;
  B.matrix[1][2] = 875421.09875;
  A.matrix[2][0] = 809.7654321;
  B.matrix[2][0] = 809.7654321;
  A.matrix[2][1] = 1.0540234543;
  B.matrix[2][1] = 1.0540234943;
  A.matrix[2][2] = -1.75;
  B.matrix[2][2] = -1.75;
  int res5 = s21_eq_matrix(&A, &B);
  ck_assert_int_eq(res5, 1);
  s21_remove_matrix(&A);
  s21_remove_matrix(&B);

  int rows8 = 2, columns8 = 3;
  s21_create_matrix(rows8, columns8, &A);
  s21_create_matrix(rows8, columns8, &B);
  A.matrix[0][0] = 1.0;
  B.matrix[0][0] = 1.0;
  A.matrix[0][1] = 0.0000007;
  B.matrix[0][1] = 0.0000007;
  A.matrix[0][2] = -16.7865444211;
  B.matrix[0][2] = -16.78654442541;
  A.matrix[1][0] = 7.08;
  B.matrix[1][0] = 7.08;
  A.matrix[1][1] = 3.11111111111;
  B.matrix[1][1] = 3.11111193;
  A.matrix[1][2] = -10.0;
  B.matrix[1][2] = -10.0;
  int res6 = s21_eq_matrix(&A, &B);
  ck_assert_int_eq(res6, 0);
  s21_remove_matrix(&A);
  s21_remove_matrix(&B);
}
END_TEST

START_TEST(test_sum_martix) {
  int rows = 2, columns = 3;
  matrix_t A = {0}, B = {0}, result = {0};
  s21_create_matrix(rows, columns, &A);
  s21_create_matrix(rows, columns, &B);
  A.matrix[0][0] = 1;
  B.matrix[0][0] = 1;
  A.matrix[0][1] = 2;
  B.matrix[0][1] = 2;
  A.matrix[0][2] = 3;
  B.matrix[0][2] = 3;
  A.matrix[1][0] = 4;
  B.matrix[1][0] = 4;
  A.matrix[1][1] = 5;
  B.matrix[1][1] = 5;
  A.matrix[1][2] = 6;
  B.matrix[1][2] = 6;
  int res = s21_sum_matrix(&A, &B, &result);
  ck_assert_int_eq(res, 0);
  ck_assert_int_eq(result.matrix[0][0], 2);
  ck_assert_int_eq(result.matrix[0][1], 4);
  ck_assert_int_eq(result.matrix[0][2], 6);
  ck_assert_int_eq(result.matrix[1][0], 8);
  ck_assert_int_eq(result.matrix[1][1], 10);
  ck_assert_int_eq(result.matrix[1][2], 12);
  s21_remove_matrix(&A);
  s21_remove_matrix(&B);
  s21_remove_matrix(&result);

  int rows5 = 2, columns5 = 2;
  s21_create_matrix(rows5, columns5, &A);
  s21_create_matrix(rows5, columns5, &B);
  A.matrix[0][0] = 1.45628;
  B.matrix[0][0] = 1.23408;
  A.matrix[0][1] = -2.43;
  B.matrix[0][1] = -1.4572;
  A.matrix[1][0] = 0.08;
  B.matrix[1][0] = 0.000123;
  A.matrix[1][1] = 5.3;
  B.matrix[1][1] = -5.1;
  int res4 = s21_sum_matrix(&A, &B, &result);
  ck_assert_int_eq(res4, 0);
  ck_assert_double_eq_tol(result.matrix[0][0], 2.69036, EPS);
  ck_assert_double_eq_tol(result.matrix[0][1], -3.8872, EPS);
  ck_assert_double_eq_tol(result.matrix[1][0], 0.080123, EPS);
  ck_assert_double_eq_tol(result.matrix[1][1], 0.2, EPS);
  s21_remove_matrix(&A);
  s21_remove_matrix(&B);
  s21_remove_matrix(&result);

  int rows2 = 2, columns2 = 3;
  int rows3 = 1, columns3 = 4;
  s21_create_matrix(rows2, columns2, &A);
  s21_create_matrix(rows3, columns3, &B);
  int res2 = s21_sum_matrix(&A, &B, &result);
  ck_assert_int_eq(res2, 2);
  s21_remove_matrix(&A);
  s21_remove_matrix(&B);

  int rows4 = 0, columns4 = -9;
  s21_create_matrix(rows4, columns4, &A);
  s21_create_matrix(rows4, columns4, &B);
  int res3 = s21_sum_matrix(&A, &B, &result);
  ck_assert_int_eq(res3, 1);
  s21_remove_matrix(&A);
  s21_remove_matrix(&B);
}
END_TEST

START_TEST(test_sub_martix) {
  int rows = 2, columns = 3;
  matrix_t A = {0}, B = {0}, result = {0};
  s21_create_matrix(rows, columns, &A);
  s21_create_matrix(rows, columns, &B);
  A.matrix[0][0] = 7;
  B.matrix[0][0] = 1;
  A.matrix[0][1] = 5;
  B.matrix[0][1] = 9;
  A.matrix[0][2] = 3;
  B.matrix[0][2] = 1;
  A.matrix[1][0] = 0;
  B.matrix[1][0] = -1;
  A.matrix[1][1] = 90;
  B.matrix[1][1] = 76;
  A.matrix[1][2] = -10;
  B.matrix[1][2] = -78;
  int res = s21_sub_matrix(&A, &B, &result);
  ck_assert_int_eq(res, 0);
  ck_assert_int_eq(result.matrix[0][0], 6);
  ck_assert_int_eq(result.matrix[0][1], -4);
  ck_assert_int_eq(result.matrix[0][2], 2);
  ck_assert_int_eq(result.matrix[1][0], 1);
  ck_assert_int_eq(result.matrix[1][1], 14);
  ck_assert_int_eq(result.matrix[1][2], 68);
  s21_remove_matrix(&A);
  s21_remove_matrix(&B);
  s21_remove_matrix(&result);

  int rows2 = 2, columns2 = 3;
  int rows3 = 1, columns3 = 4;
  s21_create_matrix(rows2, columns2, &A);
  s21_create_matrix(rows3, columns3, &B);
  int res2 = s21_sub_matrix(&A, &B, &result);
  ck_assert_int_eq(res2, 2);
  s21_remove_matrix(&A);
  s21_remove_matrix(&B);

  int rows4 = 0, columns4 = -9;
  s21_create_matrix(rows4, columns4, &A);
  s21_create_matrix(rows4, columns4, &B);
  int res3 = s21_sub_matrix(&A, &B, &result);
  ck_assert_int_eq(res3, 1);
  s21_remove_matrix(&A);
  s21_remove_matrix(&B);

  int rows5 = 2, columns5 = 2;
  s21_create_matrix(rows5, columns5, &A);
  s21_create_matrix(rows5, columns5, &B);
  A.matrix[0][0] = 1.45628;
  B.matrix[0][0] = 1.23408;
  A.matrix[0][1] = -2.43;
  B.matrix[0][1] = -1.4572;
  A.matrix[1][0] = 0.08;
  B.matrix[1][0] = 0.000123;
  A.matrix[1][1] = 5.3;
  B.matrix[1][1] = -5.1;
  int res4 = s21_sub_matrix(&A, &B, &result);
  ck_assert_int_eq(res4, 0);
  ck_assert_double_eq_tol(result.matrix[0][0], 0.2222, EPS);
  ck_assert_double_eq_tol(result.matrix[0][1], -0.9728, EPS);
  ck_assert_double_eq_tol(result.matrix[1][0], 0.079877, EPS);
  ck_assert_double_eq_tol(result.matrix[1][1], 10.4, EPS);
  s21_remove_matrix(&A);
  s21_remove_matrix(&B);
  s21_remove_matrix(&result);
}
END_TEST

START_TEST(test_mult_number) {
  int rows = 2, columns = 3, number = 2;
  matrix_t A = {0}, result = {0};
  s21_create_matrix(rows, columns, &A);
  A.matrix[0][0] = 7;
  A.matrix[0][1] = 5;
  A.matrix[0][2] = 3;
  A.matrix[1][0] = 0;
  A.matrix[1][1] = 90;
  A.matrix[1][2] = -10;
  int res = s21_mult_number(&A, number, &result);
  ck_assert_int_eq(res, 0);
  ck_assert_int_eq(result.matrix[0][0], 14);
  ck_assert_int_eq(result.matrix[0][1], 10);
  ck_assert_int_eq(result.matrix[0][2], 6);
  ck_assert_int_eq(result.matrix[1][0], 0);
  ck_assert_int_eq(result.matrix[1][1], 180);
  ck_assert_int_eq(result.matrix[1][2], -20);
  s21_remove_matrix(&A);
  s21_remove_matrix(&result);

  int rows2 = -10, columns2 = 3, number2 = -1;
  s21_create_matrix(rows2, columns2, &A);
  int res2 = s21_mult_number(&A, number2, &result);
  ck_assert_int_eq(res2, 1);
  s21_remove_matrix(&A);

  int rows3 = 10, columns3 = -3, number3 = -1;
  s21_create_matrix(rows3, columns3, &A);
  int res3 = s21_mult_number(&A, number3, &result);
  ck_assert_int_eq(res3, 1);
  s21_remove_matrix(&A);

  int rows4 = 2, columns4 = 2;
  double number4 = 0.2;
  s21_create_matrix(rows4, columns4, &A);
  A.matrix[0][0] = 7.1;
  A.matrix[0][1] = -0.000002;
  A.matrix[1][0] = 0.34;
  A.matrix[1][1] = 9.02;
  int res4 = s21_mult_number(&A, number4, &result);
  ck_assert_int_eq(res4, 0);
  ck_assert_double_eq_tol(result.matrix[0][0], 1.42, EPS);
  ck_assert_double_eq_tol(result.matrix[0][1], -0.0000004, EPS);
  ck_assert_double_eq_tol(result.matrix[1][0], 0.068, EPS);
  ck_assert_double_eq_tol(result.matrix[1][1], 1.804, EPS);
  s21_remove_matrix(&A);
  s21_remove_matrix(&result);
}
END_TEST

START_TEST(test_mult_matrix) {
  int rowsA = 2, columnsA = 3;
  int rowsB = 3, columnsB = 2;
  matrix_t A = {0}, B = {0}, result = {0};
  s21_create_matrix(rowsA, columnsA, &A);
  s21_create_matrix(rowsB, columnsB, &B);
  A.matrix[0][0] = 1;
  A.matrix[0][1] = 2;
  A.matrix[0][2] = 3;
  A.matrix[1][0] = 4;
  A.matrix[1][1] = 5;
  A.matrix[1][2] = 6;
  B.matrix[0][0] = 7;
  B.matrix[0][1] = 8;
  B.matrix[1][0] = 9;
  B.matrix[1][1] = 10;
  B.matrix[2][0] = 11;
  B.matrix[2][1] = 12;
  int res = s21_mult_matrix(&A, &B, &result);
  ck_assert_int_eq(res, 0);
  ck_assert_double_eq(result.matrix[0][0], 58.0);
  ck_assert_double_eq(result.matrix[0][1], 64.0);
  ck_assert_double_eq(result.matrix[1][0], 139.0);
  ck_assert_double_eq(result.matrix[1][1], 154.0);
  s21_remove_matrix(&A);
  s21_remove_matrix(&B);
  s21_remove_matrix(&result);

  int rowsA2 = 2, columnsA2 = 3;
  int rowsB2 = 4, columnsB2 = 2;
  s21_create_matrix(rowsA2, columnsA2, &A);
  s21_create_matrix(rowsB2, columnsB2, &B);
  int res2 = s21_mult_matrix(&A, &B, &result);
  ck_assert_int_eq(res2, 2);
  s21_remove_matrix(&A);
  s21_remove_matrix(&B);

  int rowsA3 = 2, columnsA3 = 3;
  int rowsB3 = 4, columnsB3 = -8;
  s21_create_matrix(rowsA3, columnsA3, &A);
  s21_create_matrix(rowsB3, columnsB3, &B);
  int res3 = s21_mult_matrix(&A, &B, &result);
  ck_assert_int_eq(res3, 1);
  s21_remove_matrix(&A);
  s21_remove_matrix(&B);
}
END_TEST

START_TEST(test_transpose) {
  int rows = 3, columns = 2;
  matrix_t A = {0}, result = {0};
  s21_create_matrix(rows, columns, &A);
  A.matrix[0][0] = 1;
  A.matrix[0][1] = 4;
  A.matrix[1][0] = 2;
  A.matrix[1][1] = 5;
  A.matrix[2][0] = 3;
  A.matrix[2][1] = 6;
  int res = s21_transpose(&A, &result);
  ck_assert_int_eq(res, 0);
  ck_assert_int_eq(result.matrix[0][0], 1);
  ck_assert_int_eq(result.matrix[0][1], 2);
  ck_assert_int_eq(result.matrix[0][2], 3);
  ck_assert_int_eq(result.matrix[1][0], 4);
  ck_assert_int_eq(result.matrix[1][1], 5);
  ck_assert_int_eq(result.matrix[1][2], 6);
  s21_remove_matrix(&A);
  s21_remove_matrix(&result);

  int rows2 = 3, columns2 = 3;
  s21_create_matrix(rows2, columns2, &A);
  A.matrix[0][0] = 4;
  A.matrix[0][1] = 2;
  A.matrix[0][2] = 4.873;
  A.matrix[1][0] = 1;
  A.matrix[1][1] = 54.08;
  A.matrix[1][2] = 3;
  A.matrix[2][0] = 0;
  A.matrix[2][1] = -6;
  A.matrix[2][2] = -1;
  int res2 = s21_transpose(&A, &result);
  ck_assert_int_eq(res2, 0);
  ck_assert_double_eq(result.matrix[0][0], 4);
  ck_assert_double_eq(result.matrix[0][1], 1);
  ck_assert_double_eq(result.matrix[0][2], 0);
  ck_assert_double_eq(result.matrix[1][0], 2);
  ck_assert_double_eq(result.matrix[1][1], 54.08);
  ck_assert_double_eq(result.matrix[1][2], -6);
  ck_assert_double_eq(result.matrix[2][0], 4.873);
  ck_assert_double_eq(result.matrix[2][1], 3);
  ck_assert_double_eq(result.matrix[2][2], -1);
  s21_remove_matrix(&A);
  s21_remove_matrix(&result);

  int rows3 = 2, columns3 = 0;
  if (!s21_create_matrix(rows3, columns3, &A)) {
    int res3 = s21_transpose(&A, &result);
    ck_assert_int_eq(res3, 1);
    s21_remove_matrix(&A);
  }
}
END_TEST

START_TEST(test_determinant) {
  matrix_t A = {0};
  double det;
  s21_create_matrix(2, 2, &A);
  A.matrix[0][0] = 1;
  A.matrix[0][1] = 2;
  A.matrix[1][0] = 3;
  A.matrix[1][1] = 4;
  int res = s21_determinant(&A, &det);
  ck_assert_int_eq(res, 0);
  ck_assert_double_eq_tol(det, -2.0, EPS);
  s21_remove_matrix(&A);

  s21_create_matrix(3, 3, &A);
  A.matrix[0][0] = 1;
  A.matrix[0][1] = 2;
  A.matrix[0][2] = 3;
  A.matrix[1][0] = 4;
  A.matrix[1][1] = 5;
  A.matrix[1][2] = 6;
  A.matrix[2][0] = 7;
  A.matrix[2][1] = 8;
  A.matrix[2][2] = 9;
  int res2 = s21_determinant(&A, &det);
  ck_assert_int_eq(res2, 0);
  ck_assert_double_eq_tol(det, 0.0, EPS);
  s21_remove_matrix(&A);

  s21_create_matrix(3, 3, &A);
  A.matrix[0][0] = 2;
  A.matrix[0][1] = 0;
  A.matrix[0][2] = 0;
  A.matrix[1][0] = 1;
  A.matrix[1][1] = 3;
  A.matrix[1][2] = 0;
  A.matrix[2][0] = 1;
  A.matrix[2][1] = 1;
  A.matrix[2][2] = 4;
  int res3 = s21_determinant(&A, &det);
  ck_assert_int_eq(res3, 0);
  ck_assert_double_eq_tol(det, 24.0, EPS);
  s21_remove_matrix(&A);

  s21_create_matrix(2, 2, &A);
  A.matrix[0][0] = 1;
  A.matrix[0][1] = -1;
  A.matrix[1][0] = 2;
  A.matrix[1][1] = 2;
  int res4 = s21_determinant(&A, &det);
  ck_assert_int_eq(res4, 0);
  ck_assert_double_eq_tol(det, 4.0, EPS);
  s21_remove_matrix(&A);

  s21_create_matrix(1, 1, &A);
  A.matrix[0][0] = 1.0;
  int res5 = s21_determinant(&A, &det);
  ck_assert_int_eq(res5, 0);
  ck_assert_double_eq_tol(det, 1.0, EPS);
  s21_remove_matrix(&A);

  s21_create_matrix(4, 4, &A);
  A.matrix[0][0] = 0;
  A.matrix[0][1] = 0;
  A.matrix[0][2] = 2;
  A.matrix[0][3] = 4;
  A.matrix[1][0] = 0;
  A.matrix[1][1] = -3;
  A.matrix[1][2] = -4;
  A.matrix[1][3] = 3;
  A.matrix[2][0] = 4;
  A.matrix[2][1] = 2;
  A.matrix[2][2] = 2;
  A.matrix[2][3] = 2;
  A.matrix[3][0] = 4;
  A.matrix[3][1] = 2;
  A.matrix[3][2] = -1;
  A.matrix[3][3] = -1;
  int res6 = s21_determinant(&A, &det);
  ck_assert_int_eq(res6, 0);
  ck_assert_double_eq_tol(det, 72.0, EPS);
  s21_remove_matrix(&A);

  s21_create_matrix(2, 3, &A);
  A.matrix[0][0] = 4;
  A.matrix[0][1] = 2;
  A.matrix[0][2] = 4.873;
  A.matrix[1][0] = 1;
  A.matrix[1][1] = 54.08;
  A.matrix[1][2] = 3;
  int res7 = s21_determinant(&A, &det);
  ck_assert_int_eq(res7, 2);
  s21_remove_matrix(&A);

  s21_create_matrix(2, 0, &A);
  int res8 = s21_determinant(&A, &det);
  ck_assert_int_eq(res8, 1);
  s21_remove_matrix(&A);
}
END_TEST

START_TEST(test_calc_complements) {
  matrix_t A = {0}, result = {0};
  s21_create_matrix(3, 3, &A);
  A.matrix[0][0] = 1;
  A.matrix[0][1] = 2;
  A.matrix[0][2] = 3;
  A.matrix[1][0] = 0;
  A.matrix[1][1] = 4;
  A.matrix[1][2] = 2;
  A.matrix[2][0] = 5;
  A.matrix[2][1] = 2;
  A.matrix[2][2] = 1;
  int res = s21_calc_complements(&A, &result);
  ck_assert_int_eq(res, OK);
  ck_assert_int_eq(result.matrix[0][0], 0);
  ck_assert_int_eq(result.matrix[0][1], 10);
  ck_assert_int_eq(result.matrix[0][2], -20);
  ck_assert_int_eq(result.matrix[1][0], 4);
  ck_assert_int_eq(result.matrix[1][1], -14);
  ck_assert_int_eq(result.matrix[1][2], 8);
  ck_assert_int_eq(result.matrix[2][0], -8);
  ck_assert_int_eq(result.matrix[2][1], -2);
  ck_assert_int_eq(result.matrix[2][2], 4);
  s21_remove_matrix(&A);
  s21_remove_matrix(&result);

  int res3 = s21_calc_complements(NULL, &result);
  ck_assert_int_eq(res3, 1);
  if (!res3) s21_remove_matrix(&result);

  s21_create_matrix(2, 3, &A);
  A.matrix[0][0] = 4;
  A.matrix[0][1] = 2;
  A.matrix[0][2] = 4.873;
  A.matrix[1][0] = 1;
  A.matrix[1][1] = 54.08;
  A.matrix[1][2] = 3;
  int res7 = s21_calc_complements(&A, &result);
  ck_assert_int_eq(res7, 2);
  s21_remove_matrix(&A);
}
END_TEST

START_TEST(test_inverse_matrix) {
  matrix_t A = {0}, result = {0};
  s21_create_matrix(3, 3, &A);
  A.matrix[0][0] = 2;
  A.matrix[0][1] = 5;
  A.matrix[0][2] = 7;
  A.matrix[1][0] = 6;
  A.matrix[1][1] = 3;
  A.matrix[1][2] = 4;
  A.matrix[2][0] = 5;
  A.matrix[2][1] = -2;
  A.matrix[2][2] = -3;
  int res = s21_inverse_matrix(&A, &result);
  ck_assert_int_eq(res, OK);
  ck_assert_int_eq(result.matrix[0][0], 1);
  ck_assert_int_eq(result.matrix[0][1], -1);
  ck_assert_int_eq(result.matrix[0][2], 1);
  ck_assert_int_eq(result.matrix[1][0], -38);
  ck_assert_int_eq(result.matrix[1][1], 41);
  ck_assert_int_eq(result.matrix[1][2], -34);
  ck_assert_int_eq(result.matrix[2][0], 27);
  ck_assert_int_eq(result.matrix[2][1], -29);
  ck_assert_int_eq(result.matrix[2][2], 24);
  s21_remove_matrix(&A);
  s21_remove_matrix(&result);

  s21_create_matrix(3, 3, &A);
  A.matrix[0][0] = -4.0;
  A.matrix[0][1] = 6.0;
  A.matrix[0][2] = 10.0;
  A.matrix[1][0] = 1.0;
  A.matrix[1][1] = 0.0;
  A.matrix[1][2] = -2.0;
  A.matrix[2][0] = 172.0;
  A.matrix[2][1] = 23.0;
  A.matrix[2][2] = 98.0;
  int res2 = s21_inverse_matrix(&A, &result);
  ck_assert_int_eq(res2, 0);
  ck_assert_double_eq_tol(result.matrix[0][0], -0.0176516, EPS);
  ck_assert_double_eq_tol(result.matrix[0][1], 0.1373753, EPS);
  ck_assert_double_eq_tol(result.matrix[0][2], 0.0046048, EPS);
  ck_assert_double_eq_tol(result.matrix[1][0], 0.1696086, EPS);
  ck_assert_double_eq_tol(result.matrix[1][1], 0.8104375, EPS);
  ck_assert_double_eq_tol(result.matrix[1][2], -0.0007675, EPS);
  ck_assert_double_eq_tol(result.matrix[2][0], -0.0088258, EPS);
  ck_assert_double_eq_tol(result.matrix[2][1], -0.4313124, EPS);
  ck_assert_double_eq_tol(result.matrix[2][2], 0.00230238, EPS);
  s21_remove_matrix(&A);
  s21_remove_matrix(&result);

  int res3 = s21_inverse_matrix(NULL, &result);
  ck_assert_int_eq(res3, 1);
  if (!res3) s21_remove_matrix(&result);

  s21_create_matrix(2, 3, &A);
  A.matrix[0][0] = 4;
  A.matrix[0][1] = 2;
  A.matrix[0][2] = 4.873;
  A.matrix[1][0] = 1;
  A.matrix[1][1] = 54.08;
  A.matrix[1][2] = 3;
  int res7 = s21_inverse_matrix(&A, &result);
  ck_assert_int_eq(res7, 2);
  s21_remove_matrix(&A);

  s21_create_matrix(1, 1, &A);
  A.matrix[0][0] = 4;
  int res8 = s21_inverse_matrix(&A, &result);
  ck_assert_int_eq(res8, 0);
  s21_remove_matrix(&A);
  s21_remove_matrix(&result);
}
END_TEST

int main(void) {
  Suite *suite = suite_create(SUITE_NAME);
  TCase *tcase = tcase_create(TCASE_NAME);
  SRunner *srunner = srunner_create(suite);
  int counted_errors = 0;

  suite_add_tcase(suite, tcase);
  tcase_add_test(tcase, test_create_matrix);
  tcase_add_test(tcase, test_remove_matrix);
  tcase_add_test(tcase, test_eq_martix);
  tcase_add_test(tcase, test_sum_martix);
  tcase_add_test(tcase, test_sub_martix);
  tcase_add_test(tcase, test_mult_number);
  tcase_add_test(tcase, test_mult_matrix);
  tcase_add_test(tcase, test_transpose);
  tcase_add_test(tcase, test_determinant);
  tcase_add_test(tcase, test_calc_complements);
  tcase_add_test(tcase, test_inverse_matrix);

  srunner_set_fork_status(srunner, CK_NOFORK);
  srunner_run_all(srunner, CK_VERBOSE);
  counted_errors = srunner_ntests_failed(srunner);
  srunner_free(srunner);

  return counted_errors == 0 ? EXIT_SUCCESS : EXIT_FAILURE;
}
