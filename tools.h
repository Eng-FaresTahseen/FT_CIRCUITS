#ifndef TOOLS_H
#define TOOLS_H
#include <stdio.h>
#include <stdlib.h>
#include <math.h>

typedef struct
{
    char elementType;
    int node1, node2;
    double value;
    double current;
    double voltage;
    double power;
} element;

typedef struct
{
    element *components;
    int size_components;
} circuit;

circuit read_circuit();
void analyse_circuit_I_Only(circuit *c);
void analyse_circuit(circuit *c);
double **find_conductance_matrix(element *components, int size_components, int n_nodes);
double *find_current_vector(element *components, int size_components, int n_nodes);
int max_node(element *c, int size);
double *solve_by_GaussElimination(double **G, int n, double *I);
void print_circuit_analysis(circuit c);
void swap_rows(int row1, int row2, double **G, double *I, int n);
void fix_swap_rows(double **G, double *I, int n);
double **get_transpose(double **matrix, int n_rows, int n_cols);
double **multiply_matrices(double **matrix1, double **matrix2, int n_rows1, int n_cols1, int n_rows2, int n_cols2);
double **add_matrices(double **matrix1, double **matrix2, int n_rows, int n_cols);
int get_number_voltage_sources(element *components, int size_components);
double **get_incidence_matrix(element *components, int n_nodes_ex0, int size_components, int n_v_sources);
short get_valid_choice(short lower, short upper);
void Analyze_circuit();
void show_manual_guide();
void Analyze_circuit_file();
void show_main_menu();
void print_resistor(int node1, int node2, double value);
void print_voltage_source(int node1, int node2, double value);
void print_current_source(int node1, int node2, double value);
void print_detailed_circuit_analysis(circuit c);
#endif