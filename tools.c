#include "tools.h"

short get_valid_choice(short lower, short upper)
{
    short choice;
    scanf("%hd", &choice);
    while ((choice > upper) || (choice < lower))
    {
        printf("\nEnter a valid choice : ");
        scanf("%hd", &choice);
    }
    return choice;
}

void show_main_menu()
{
    printf("\n________________________________________________________\n");
    printf("\n                      FT Circuits                       ");
    printf("\n                        Analyzer                        \n");
    printf("\n________________________________________________________\n");
    printf("[1] Manual Guide\n");
    printf("[2] Analyze circuit\n");
    printf("[3] Analyze circuit file\n");
    printf("[4] Exit");
    printf("\n________________________________________________________\n");

    short choice = get_valid_choice(1, 4);

    switch (choice)
    {
    case 1:
        show_manual_guide();
        break;
    case 2:
        Analyze_circuit();
        break;
    case 3:
        Analyze_circuit_file();
        break;
    case 4:
        printf("\n---------------------------------------------\n");
        printf("        The program has ended successfully\n");
        printf("---------------------------------------------\n");
        break;
    default:
        break;
    }
}

void Analyze_circuit()
{
    circuit c = read_circuit();
    analyse_circuit(&c);
    printf("\n---------------------------------------------\n");
    printf("The circuit has been analyzed successfully!");
    printf("\n---------------------------------------------\n");
    short choice_detail;
    printf("Do you want to see the detailed analysis? \n[1] Yes\n[2] No\n");
    choice_detail = get_valid_choice(1, 2);
    if (choice_detail == 1)
    {
        print_detailed_circuit_analysis(c);
    }
    else
    {
        print_circuit_analysis(c);
    }

    free(c.components);

    printf("\n---------------------------------------------\n");

    printf("\n[1] Go back to the Main Menu\n");
    printf("[2] Analyze another circuit\n");
    printf("[3] Analyze circuit file \n");
    printf("[4] Exit the program\n");
    short choice = get_valid_choice(1, 4);
    switch (choice)
    {
    case 1:
        show_main_menu();
        break;
    case 2:
        Analyze_circuit();
        break;
    case 3:
        Analyze_circuit_file();
        break;
    case 4:
        printf("\n---------------------------------------------\n");
        printf("        The program has ended successfully\n");
        printf("---------------------------------------------\n");
        break;
    default:
        break;
    }
}

void Analyze_circuit_file()
{
    circuit c;
    c.components = NULL;
    c.size_components = 0;

    char file_name[100], full_path[200];
    printf("\nNote : The file must be inside the 'circuits_files' folder");
    printf("\nEnter the file name : ");
    scanf("%s", file_name);
    sprintf(full_path, "circuits_files/%s", file_name);
    FILE *file = fopen(full_path, "r");
    if (!file)
    {
        printf("Error has occurred during opening the file ..!\n");
        return;
    }
    char line[100];
    int n = 0;
    while (fgets(line, sizeof(line), file))
    {
        if (line[0] == '\n')
            continue;
        n++;
    }
    rewind(file);
    c.size_components = n;

    // validate the number of elements

    if (n <= 1)
    {
        printf("Error: Number of elements is not valid.\n");
        fclose(file);
        return;
    }

    c.components = malloc(n * sizeof(element));
    int i = 0;
    while (fgets(line, sizeof(line), file))
    {
        if (line[0] == '\n')
            continue;
        sscanf(line, " %c %d %d %lf",
               &c.components[i].elementType,
               &c.components[i].node1,
               &c.components[i].node2,
               &c.components[i].value);
        i++;
    }
    fclose(file);

    // validate if there is a power source
    int has_power_source = 0;
    for (int i = 0; i < n; i++)
    {
        if (c.components[i].elementType == 'V' || c.components[i].elementType == 'I')
        {
            has_power_source = 1;
            break;
        }
    }
    if (!has_power_source)
    {
        printf("Error: Circuit must contain a power source (V or I).\n");
        free(c.components);
        return;
    }

    analyse_circuit(&c);
    printf("\n---------------------------------------------\n");
    printf("The circuit has been analyzed successfully!");
    printf("\n---------------------------------------------\n");
    short choice_detail;
    printf("Do you want to see the detailed analysis? \n[1] Yes\n[2] No\n");
    choice_detail = get_valid_choice(1, 2);
    if (choice_detail == 1)
    {
        print_detailed_circuit_analysis(c);
    }
    else
    {
        print_circuit_analysis(c);
    }
    free(c.components);

    printf("\n---------------------------------------------\n");

    printf("\n[1] Go back to the Main Menu\n");
    printf("[2] Analyze circuit\n");
    printf("[3] Analyze another circuit file \n");
    printf("[4] Exit the program\n");
    short choice = get_valid_choice(1, 4);
    switch (choice)
    {
    case 1:
        show_main_menu();
        break;
    case 2:
        Analyze_circuit();
        break;
    case 3:
        Analyze_circuit_file();
        break;
    case 4:
        printf("\n---------------------------------------------\n");
        printf("        The program has ended successfully\n");
        printf("---------------------------------------------\n");
        break;
    default:
        break;
    }
}

void show_manual_guide()
{
    printf("=============================================\n");
    printf("        FT CIRCUIT ANALYZER USER GUIDE\n");
    printf("=============================================\n\n");

    printf("1) NODE NUMBERING RULES:\n");
    printf("---------------------------------------------\n");
    printf("- Node 0 is ALWAYS the ground (reference).\n");
    printf("- All other nodes must be numbered starting from 1.\n");
    printf("- Node numbers must be non-negative integers.\n");
    printf("- Example: 0, 1, 2, 3, ...\n\n");

    printf("2) ELEMENT INPUT FORMAT:\n");
    printf("---------------------------------------------\n");
    printf("Each element must be entered in the form:\n\n");
    printf("  <Type> <Node1> <Node2> <Value>\n\n");

    printf("Where:\n");
    printf("- <Type>  : Element type (R, I, V)\n");
    printf("- <Node1> : First connection node\n");
    printf("- <Node2> : Second connection node\n");
    printf("- <Value> : Element value\n\n");

    printf("3) SUPPORTED ELEMENT TYPES:\n");
    printf("---------------------------------------------\n");
    printf("R : Resistor (Ohms)\n");
    printf("I : Current Source (Amperes)\n");
    printf("V : Voltage Source (Volts)\n\n");

    printf("4) SIGN CONVENTIONS:\n");
    printf("---------------------------------------------\n");
    printf("- Resistor R between Node1 and Node2 has no polarity.\n");
    printf("- Current source I flows FROM Node1 TO Node2.\n");
    printf("- Voltage source V has positive terminal at Node2\n");
    printf("  and negative terminal at Node1.\n\n");

    printf("5) EXAMPLES:\n");
    printf("---------------------------------------------\n");
    printf("R 1 2 100     -> Resistor of 100 Ohms between nodes 1 and 2\n");
    printf("I 1 0 2       -> 2A current source from node 1 to ground\n");
    printf("V 2 0 5       -> 5V voltage source (node ground positive, 2 negative)\n\n");

    printf("6) IMPORTANT NOTES:\n");
    printf("---------------------------------------------\n");
    printf("- Do NOT skip node numbers.\n");
    printf("- Ground (node 0) must be connected to the circuit.\n");
    printf("- All values must be real numbers.\n");
    printf("- The solver uses Modified Nodal Analysis (MNA).\n\n");

    printf("=============================================\n");
    printf("          END OF USER GUIDE\n");
    printf("=============================================\n");

    printf("[1] Go back to the Main Menu\n");
    printf("[2] Exit the program\n");

    short choice = get_valid_choice(1, 2);
    switch (choice)
    {
    case 1:
        show_main_menu();
        break;
    case 2:
        printf("\n---------------------------------------------\n");
        printf("        The program has ended successfully\n");
        printf("---------------------------------------------\n");
        break;
    default:
        break;
    }
}

circuit read_circuit()
{
    printf("how many elements exist in the circuit : ");
    int n;
    scanf("%d", &n);
    if (n <= 1)
    {
        printf("Error: Number of elements is not valid.\n");
        exit(EXIT_FAILURE);
    }
    element *components = malloc(n * sizeof(element));
    printf("Enter the elements in the separated form with nodes given that node 0 is the ground \nfor example : R 1 2 200 -> that respresents a resistor of 200 ohms between nodes 1 and 2 :\n");
    for (int i = 0; i < n; i++)
    {
        scanf(" %c %d %d %lf", &(components[i].elementType), &(components[i].node1), &(components[i].node2), &(components[i].value));
    }

    // validate if there is a power source
    int has_power_source = 0;
    for (int i = 0; i < n; i++)
    {
        if (components[i].elementType == 'V' || components[i].elementType == 'I')
        {
            has_power_source = 1;
            break;
        }
    }
    if (!has_power_source)
    {
        printf("Error: Circuit must contain a power source (V or I).\n");
        exit(EXIT_FAILURE);
    }

    circuit c;
    c.components = components;
    c.size_components = n;
    return c;
}

int max_node(element *c, int size)
{
    int max = 0;
    for (int i = 0; i < size; i++)
    {
        if (c[i].node1 > max)
            max = c[i].node1;
        if (c[i].node2 > max)
            max = c[i].node2;
    }
    return max;
}

double **find_conductance_matrix(element *components, int size_components, int n_nodes)
{
    double **G = malloc((n_nodes - 1) * sizeof(double *));
    for (int i = 0; i < n_nodes - 1; i++)
    {
        G[i] = calloc((n_nodes - 1), sizeof(double));
    }
    for (int i = 0; i < size_components; i++)
    {
        if (components[i].elementType != 'R')
            continue;
        int x = components[i].node1, y = components[i].node2;
        double current_cond = (1.00 / components[i].value);
        if (x > 0)
            G[x - 1][x - 1] += current_cond;
        if (y > 0)
            G[y - 1][y - 1] += current_cond;
        if ((x > 0) && (y > 0))
        {
            G[x - 1][y - 1] -= current_cond;
            G[y - 1][x - 1] -= current_cond;
        }
    }
    return G;
}

double *find_current_vector(element *components, int size_components, int n_nodes)
{
    double *I = calloc(n_nodes - 1, sizeof(double));
    for (int i = 0; i < size_components; i++)
    {
        if (components[i].elementType != 'I')
            continue;
        int from_node = components[i].node1 - 1, to_node = components[i].node2 - 1;
        double current_value = components[i].value;
        if (from_node >= 0)
            I[from_node] -= current_value;
        if (to_node >= 0)
            I[to_node] += current_value;
    }
    return I;
}

double *solve_by_GaussElimination(double **G, int n, double *I)
{
    double *V = calloc(n, sizeof(double));
    fix_swap_rows(G, I, n);

    // Forward Elimination
    for (int row = 0; row < n; row++)
    {
        for (int col = 0; col < n; col++)
        {
            if (G[row][col] == 0)
                continue;
            for (int i = row + 1; i < n; i++)
            {
                double factor = -(G[i][col] / G[row][col]);
                for (int k = col; k < n; k++)
                {
                    G[i][k] += G[row][k] * factor;
                }
                I[i] += factor * I[row];
            }
            break;
        }
    }

    // Backward substitution
    for (int row = n - 1; row >= 0; row--)
    {
        double sum = I[row];
        for (int col = row + 1; col < n; col++)
        {
            sum -= G[row][col] * V[col];
        }
        V[row] = sum / G[row][row];
    }
    return V;
}

void print_resistor(int node1, int node2, double value)
{
    printf("          %.2f Ω   \n", value);
    printf("Node %d o---^^^---o Node %d\n", node1, node2);
}

void print_voltage_source(int node1, int node2, double value)
{
    printf("           %.2f V   \n", value);
    printf("Node %d o---[- +]---o Node %d\n", node1, node2);
}

void print_current_source(int node1, int node2, double value)
{
    printf("Node %d o-->-(I %.2f A)->--o Node %d\n", node1, value, node2);
}

void print_detailed_circuit_analysis(circuit c)
{
    printf("Detailed Circuit Analysis:\n");
    for (int i = 0; i < c.size_components; i++)
    {
        switch (c.components[i].elementType)
        {
        case 'R':
            print_resistor(c.components[i].node1, c.components[i].node2, c.components[i].value);
            break;
        case 'V':
            print_voltage_source(c.components[i].node1, c.components[i].node2, c.components[i].value);
            break;
        case 'I':
            print_current_source(c.components[i].node1, c.components[i].node2, c.components[i].value);
            break;
        default:
            break;
        }
        printf("  Voltage: %.5lf V, Current: %.5lf A, Power: %.5lf W\n\n",
               c.components[i].voltage,
               c.components[i].current,
               c.components[i].power);
    }
}

void print_circuit_analysis(circuit c)
{
    printf("Element Analysis:\n");
    for (int i = 0; i < c.size_components; i++)
    {
        printf("Element %d: Type: %c, Nodes: (%d, %d), Value: %.3lf, Voltage: %.3lf V, Current: %.3lf A, Power: %.3lf W\n",
               i + 1,
               c.components[i].elementType,
               c.components[i].node1,
               c.components[i].node2,
               c.components[i].value,
               c.components[i].voltage,
               c.components[i].current,
               c.components[i].power);
    }
}

void swap_rows(int row1, int row2, double **G, double *I, int n)
{
    double temp;
    for (int i = 0; i < n; i++)
    {
        temp = G[row1][i];
        G[row1][i] = G[row2][i];
        G[row2][i] = temp;
    }
    temp = I[row1];
    I[row1] = I[row2];
    I[row2] = temp;
}

void fix_swap_rows(double **G, double *I, int n)
{
    for (int i = 0; i < n; i++)
    {
        int k = i + 1;
        if (G[i][i] == 0)
        {
            while (k < n)
            {
                if (G[k][i] != 0)
                {
                    swap_rows(k, i, G, I, n);
                    break;
                }
                k++;
            }
        }
    }
}

double **get_transpose(double **matrix, int n_rows, int n_cols)
{
    double **transpose = malloc(n_cols * sizeof(double *));
    for (int i = 0; i < n_cols; i++)
    {
        transpose[i] = malloc(n_rows * sizeof(double));
    }
    for (int i = 0; i < n_rows; i++)
    {
        for (int j = 0; j < n_cols; j++)
        {
            transpose[j][i] = matrix[i][j];
        }
    }
    return transpose;
}

double **add_matrices(double **matrix1, double **matrix2, int n_rows, int n_cols)
{
    double **added = malloc(n_rows * sizeof(double *));
    for (int i = 0; i < n_rows; i++)
    {
        added[i] = malloc(n_cols * sizeof(double));
    }
    for (int i = 0; i < n_rows; i++)
    {
        for (int j = 0; j < n_cols; j++)
        {
            added[i][j] = matrix1[i][j] + matrix2[i][j];
        }
    }
    return added;
}

double **multiply_matrices(double **matrix1, double **matrix2, int n_rows1, int n_cols1, int n_rows2, int n_cols2)
{
    if (n_cols1 != n_rows2)
    {
        printf("\nError has occurred ..!\n");
        return NULL;
    }

    double **res = malloc(n_rows1 * sizeof(double *));
    for (int i = 0; i < n_rows1; i++)
    {
        res[i] = malloc(n_cols2 * sizeof(double));
    }

    for (int row = 0; row < n_rows1; row++)
    {
        for (int col = 0; col < n_cols2; col++)
        {
            res[row][col] = 0;
            for (int k = 0; k < n_rows2; k++)
            {
                res[row][col] += matrix1[row][k] * matrix2[k][col];
            }
        }
    }
    return res;
}

int get_number_voltage_sources(element *components, int size_components)
{
    int m = 0;
    for (int i = 0; i < size_components; i++)
    {
        if (components[i].elementType == 'V')
        {
            m++;
        }
    }
    return m;
}

double **get_incidence_matrix(element *components, int n_nodes_ex0, int size_components, int n_v_sources)
{
    double **B = malloc(n_nodes_ex0 * sizeof(double *));
    for (int i = 0; i < n_nodes_ex0; i++)
    {
        B[i] = calloc(n_v_sources, sizeof(double));
    }
    int j = 0;
    for (int i = 0; i < size_components; i++)
    {
        if (components[i].elementType != 'V')
            continue;
        int n1 = components[i].node1 - 1;
        int n2 = components[i].node2 - 1;

        if (n1 >= 0)
            B[n1][j] = -1;
        if (n2 >= 0)
            B[n2][j] = 1;

        j++;
    }
    return B;
}

double *get_voltage_source_values(circuit *c, int m)
{
    double *E = malloc(m * sizeof(double));
    int k = 0;
    for (int i = 0; i < c->size_components; i++)
    {
        if (c->components[i].elementType == 'V')
        {
            E[k++] = c->components[i].value;
        }
    }
    return E;
}

void analyse_circuit_I_Only(circuit *c)
{
    int n_nodes = max_node(c->components, c->size_components) + 1;
    double v_nodes[n_nodes], *temp, **G = find_conductance_matrix(c->components, c->size_components, n_nodes),
                                    *I = find_current_vector(c->components, c->size_components, n_nodes);
    // validate if G is singular
    for (int i = 0; i < n_nodes - 1; i++)
    {
        if (G[i][i] == 0)
        {
            printf("\nThe circuit is invalid (The circuit is open) ..!\n");
            exit(EXIT_FAILURE);
        }
    }

    v_nodes[0] = 0;
    temp = solve_by_GaussElimination(G, n_nodes - 1, I);
    for (int i = 0; i < n_nodes - 1; i++)
    {
        v_nodes[i + 1] = temp[i];
    }
    for (int i = 0; i < c->size_components; i++)
    {
        c->components[i].voltage = v_nodes[c->components[i].node2] - v_nodes[c->components[i].node1];
        if (c->components[i].elementType == 'R')
        {
            c->components[i].current = (c->components[i].voltage) / (c->components[i].value);
        }
        else if (c->components[i].elementType == 'I')
        {
            c->components[i].current = c->components[i].value;
        }
        c->components[i].power = (c->components[i].current) * (c->components[i].voltage);
    }

    // free matrices
    for (int i = 0; i < n_nodes - 1; i++)
    {
        free(G[i]);
    }
    free(G);
    free(I);
}

void analyse_circuit_I_V(circuit *c, int m)
{
    int n = max_node(c->components, c->size_components);
    double **M = malloc((n + m) * sizeof(double *));
    for (int i = 0; i < n + m; i++)
    {
        M[i] = calloc(n + m, sizeof(double));
    }
    double **G = find_conductance_matrix(c->components, c->size_components, n + 1),
           *I = find_current_vector(c->components, c->size_components, n + 1),
           **B = get_incidence_matrix(c->components, n, c->size_components, m),
           **C = get_transpose(B, n, m),
           *E = get_voltage_source_values(c, m);

    // validate if G is singular
    for (int i = 0; i < n; i++)
    {
        if (G[i][i] == 0)
        {
            printf("\nThe circuit is invalid (The circuit is open) ..!\n");
            exit(EXIT_FAILURE);
        }
    }

    // M * V = F

    // construct M
    for (int i = 0; i < n; i++)
    {
        for (int j = 0; j < n; j++)
        {
            M[i][j] = G[i][j];
        }
        for (int k = n; k < n + m; k++)
        {
            M[i][k] = B[i][k - n];
        }
    }
    for (int i = n; i < n + m; i++)
    {
        for (int j = 0; j < n; j++)
        {
            M[i][j] = C[i - n][j];
        }
    }

    // construct F
    double F[n + m];
    for (int i = 0; i < n; i++)
    {
        F[i] = I[i];
    }
    for (int j = 0; j < m; j++)
    {
        F[j + n] = E[j];
    }

    // Solve for the V
    double *V = solve_by_GaussElimination(M, n + m, F);

    double v_nodes[n + 1];
    v_nodes[0] = 0;
    for (int i = 0; i < n; i++)
    {
        v_nodes[i + 1] = V[i];
    }
    int t = n;
    for (int i = 0; i < c->size_components; i++)
    {
        if (c->components[i].elementType == 'R')
        {
            c->components[i].voltage = v_nodes[c->components[i].node2] - v_nodes[c->components[i].node1];
            c->components[i].current = (c->components[i].voltage) / (c->components[i].value);
        }
        else if (c->components[i].elementType == 'I')
        {
            c->components[i].voltage = v_nodes[c->components[i].node2] - v_nodes[c->components[i].node1];
            c->components[i].current = c->components[i].value;
        }
        else if (c->components[i].elementType == 'V')
        {
            c->components[i].voltage = c->components[i].value;
            c->components[i].current = V[t++];
        }
        c->components[i].power = (c->components[i].current) * (c->components[i].voltage);
    }

    // free the matrix

    for (int i = 0; i < n + m; i++)
    {
        free(M[i]);
        if (i < n)
        {
            free(G[i]);
            free(B[i]);
        }
        if (i >= n)
        {
            free(C[i - n]);
        }
    }
    free(M);
    free(G);
    free(B);
    free(C);
    free(I);
    free(E);
}

void analyse_circuit(circuit *c)
{
    int m = get_number_voltage_sources(c->components, c->size_components);
    if (m == 0)
        analyse_circuit_I_Only(c);
    else
        analyse_circuit_I_V(c, m);
}