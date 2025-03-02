#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include <string>
#include <cmath>
#include <iterator>
#include <Eigen/Eigenvalues>
#include "matplotlibcpp.h"
#include <cmath>

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

using namespace std;
using namespace Eigen;
namespace plt = matplotlibcpp;

//Metodos para el programador: Impresion de matrices y obtencion de columnas
void imprimirMatriz(const vector<vector<double>>& matriz) {
    int row = matriz.size();
    int col = matriz[0].size();


    for (int i = 0; i < row; i++) {
        for (int j = 0; j < col; j++) {
            cout << matriz[i][j] << "\t";
        }
        cout << endl;
    }
}
void imprimirMatriz(const vector<vector<string>>& matriz) {
    int row = matriz.size();
    int col = matriz[0].size();


    for (int i = 0; i < row; i++) {
        for (int j = 0; j < col; j++) {
            cout << matriz[i][j] << "\t";
        }
        cout << endl;
    }
}
pair<vector<double>, vector<double>> getColumnas(const vector<vector<double>>& matriz, int columna1, int columna2) {
    vector<double> columna1V,columna2V;
	for (int i = 0; i < matriz.size(); i++) {
		columna1V.push_back(matriz[i][columna1]);
		columna2V.push_back(matriz[i][columna2]);
	}
	return { columna1V,columna2V };
}

//Lectura de archivos CSV
vector<vector<double>> leerCSV(const string& nombreArchivo,vector<string>& y_values, vector<string>& x_values) {
    ifstream archivo(nombreArchivo);
    vector<vector<double>> datos;
    string linea;

    if (!archivo.is_open()) {
        cerr << "Error al abrir el archivo" << endl;
        return datos;
    }

    //XValues
    if (getline(archivo, linea)) {
        stringstream ss(linea);
        string materia;
        getline(ss, materia, ';'); // Saltar primer valor vacío
        while (getline(ss, materia, ';')) {
            x_values.push_back(materia);
        }
    }

    while (getline(archivo, linea)) {
        stringstream ss(linea);
        string valor;
        vector<double> fila;

        //YValues
        if (getline(ss, valor, ';')) {
            y_values.push_back(valor);
        }

        while (getline(ss, valor, ';')) {
            // cambia ',' por '.' para decimales
            for (char& c : valor) {
                if (c == ',') c = '.';
            }
            fila.push_back(stod(valor));
        }
        datos.push_back(fila);
    }

    archivo.close();
    return datos;
}

//Metodo para centrar y reducir la tabla original de datos X (Paso 1)
vector<vector<double>> estandarizar(const vector<vector<double>>& matriz) {
    vector<double> media;
    vector<double> desviacionEstandar;
    vector<vector<double>> matrizE;

    int row = matriz.size();
    int col = matriz[0].size();

    double sum = 0;
    double avg = 0;

    //media
    for (int i = 0; i < col; i++) {
        for (int j = 0; j < row; j++) {
            sum += matriz[j][i];
        }
        avg = sum / row;
        media.push_back(avg);
        sum = 0;
    }

    double sum2 = 0;

    //desviacion estandar
    for (int i = 0; i < col; i++) {
        for (int j = 0; j < row; j++) {
            sum2 += pow(matriz[j][i] - media[i],2);
        }
        avg = sqrt(sum2 / (row ));
		//cout << avg << endl;

        desviacionEstandar.push_back(avg);
        sum2 = 0;
    }

    double x = 0;
    vector<double> temp;
    //estandarizar
    for (int i = 0; i < row; i++) {
        for (int j = 0; j < col; j++) {
            if(desviacionEstandar[j] != 0){
                x = (matriz[i][j] - media[j]) / desviacionEstandar[j];
                temp.push_back(x);
            }
            else {
                temp.push_back(0);
            }
            
        }
        matrizE.push_back(temp);
        temp.clear();
    }

    return matrizE;
}

//Metodo para calcular la matriz de correlaciones (Paso 2)
vector<vector<double>> correlaciones(const vector<vector<double>>& matriz) {
    vector<vector<double>> matrizCorrelacion;

    vector<double> temp;

    int row = matriz.size();
    int col = matriz[0].size();

    double sum = 0;
    double val = 0;

    for (int i = 0; i < col; i++) {
        for (int j = 0; j < col; j++) {
            for (int k = 0; k < row; k++) {
                sum += matriz[k][i] * matriz[k][j];
            }
            val = sum / (row );
            temp.push_back(val);
            sum = 0;
        }
        matrizCorrelacion.push_back(temp);
        temp.clear();
    }

    return matrizCorrelacion;
}

//Metodo para calcular los valores propios y vectores propios y ordenarlos de mayor a menor (Paso 3 y 4)
pair<vector<double>, vector<vector<double>>> calcularValoresVectoresPropios(const vector<vector<double>>& matriz) {
    int row = matriz.size();  

    // convertir de tipo vector a MatrixXd
    MatrixXd mat(row, row);
    for (int i = 0; i < row; i++) {
        for (int j = 0; j < row; j++) {
            mat(i, j) = matriz[i][j];
        }
    }

	SelfAdjointEigenSolver<MatrixXd> solver(mat); //Clase para calcular valores propios y vectores propios

	VectorXd valoresPropios = solver.eigenvalues(); //valores propios  
	MatrixXd vectoresPropios = solver.eigenvectors();//vectores propios
	//cout << "Eigenvectors: " << solver.eigenvectors() << endl;

    VectorXd valoresOrdenados = valoresPropios.reverse();
    vector<double> valoresOV;
    for (int i = 0; i < row; i++) {
        valoresOV.push_back(valoresOrdenados(i));
    }

    vector<vector<double>> vectoresVector(row, vector<double>(row));

    for (int i = 0; i < row; i++) {
        for (int j = 0; j < row; j++) {
            vectoresVector[i][j] = vectoresPropios(i,j);
        }
    }
    
    //ordenar vectores propios de mayor a menor
    vector<vector<double>> vectoresOrdenados;
    vector<double> temp;

    for (int i = 0; i < row; i++) {
        for (int j = row-1; j >= 0; j--) {
            temp.push_back(vectoresVector[i][j]);
        }
        vectoresOrdenados.push_back(temp);
        temp.clear();
    }

    return{valoresOV, vectoresOrdenados};
}

//Metodo para calcular las componentes principales (Paso 5)
vector<vector<double>> componentes(const vector<vector<double>>& matrizX, const vector<vector<double>>& matrizV) {
    vector<vector<double>> matrizC;
    vector<double> temp;

    int row = matrizX.size();
    int col = matrizV[0].size();
	int com_value = matrizX[0].size();//valor comun de matrices
	double calculo = 0;//variable para los calculos de cada elemento de la matriz

    for (int i = 0; i < row; i++)
    {
        for (int j = 0; j < col; j++)
        {
            for (int k = 0; k < com_value; k++)
            {
				calculo += matrizX[i][k] * matrizV[k][j];
            }
			temp.push_back(calculo);
            calculo = 0;

        }
		matrizC.push_back(temp);
		temp.clear();
    }

    return matrizC;
}

//Metodo para calcular las calidades de los individuos (Paso 6)
vector<vector<double>> calidades(const vector<vector<double>>& matrizX, const vector<vector<double>>& matrizC) {
    vector<vector<double>> matrizQ;
    vector<double> temp;

    int row = matrizX.size();
    int col = matrizX[0].size();
	double calculo = 0, suma = 0;//variable para los calculos de cada elemento de la matriz y suma de elementos 

    for (int i = 0; i < row; i++)
    {
        //Sumatoria
        suma = 0;
        for (int k = 0; k < col; k++)
        {
            suma += (matrizX[i][k]) * (matrizX[i][k]);
        }

        for (int j = 0; j < col; j++)
        {
			calculo = (matrizC[i][j] * matrizC[i][j]) / suma;
			temp.push_back(calculo);
        }
		matrizQ.push_back(temp);
		temp.clear();
    }
    return matrizQ;
}

//Metodo para calcular las coordenadas de las variables (Paso 7)
vector<vector<double>> coordenadasVariables(const vector<double>& valoresPropios, const vector<vector<double>>& vectoresPropios) {
    int m = valoresPropios.size();
    vector<vector<double>> matrizT(m, vector<double>(m));
    //Matriz de coordenadas: T = V / sqrt(lambda)
    for (int i = 0; i < m; i++) {
        double lambda = valoresPropios[i];
        for (int j = 0; j < m; j++) {
            matrizT[i][j] = vectoresPropios[i][j] * sqrt(lambda);
        }
    }
    return matrizT;
}

//Metodo para calcular las calidades de las variables (Paso 8)
vector<vector<double>> calidadesVariables(const vector<vector<double>>& matrizT, const vector<double>& valoresPropios) {
    int filas = matrizT.size();
    int columnas = matrizT[0].size();

    vector<vector<double>> matrizS(filas, vector<double>(columnas));

    for (int i = 0; i < filas; i++) {
        for (int j = 0; j < columnas; j++) {
            matrizS[i][j] = (matrizT[i][j] * matrizT[i][j]) / valoresPropios[j];
        }
    }
    return matrizS;
}

//Metodo para calcular el vector de inercias (Paso 9)
vector<double> vectorInercias(const vector<double>& valoresPropios) {
    double sumaValores = 0;
    for (double lambda : valoresPropios) {
        sumaValores += lambda;
    }

    vector<double> inercias;
    for (double lambda : valoresPropios) {
        inercias.push_back(lambda * 100 / sumaValores);
    }
    return inercias;
}

//Metodos para graficar
void PlanoPrincipal(vector<string>& etiquetas, vector<vector<double>> componentesP, int columna1, int columna2) {

    pair<vector<double>, vector<double>> columnas = getColumnas(componentesP, columna1, columna2);
    vector<double> C1 = columnas.first;
    vector<double> C2 = columnas.second;

	plt::scatter(C1, C2, 100);//Grafico de dispersion

    for (size_t i = 0; i < C1.size(); ++i) {
        plt::annotate(etiquetas[i], C1[i], C2[i]);
    }

    plt::xlabel(("Componente "+ to_string(columna1)));
    plt::ylabel(("Componente "+ to_string(columna2)));
    plt::title("Plano Principal");
    plt::grid(true);//Los cuadriculas del grafico
    plt::show();
}
void Circulo_Correlacion(vector<string>& etiquetas, vector<vector<double>> matrizT, int columna1, int columna2, vector<double>vectorI) {

    pair<vector<double>, vector<double>> columnas = getColumnas(matrizT, columna1, columna2);
    vector<double> T1 = columnas.first;
    vector<double> T2 = columnas.second;
    double x = vectorI[columna1], y = vectorI[columna2];

    plt::figure_size(600, 600);

    plt::plot({ -1.0, 1.0 }, { 0, 0 }, "k--");
    plt::plot({ 0, 0 }, { -1.0, 1.0 }, "k--");

    //Grafico de las flechas
	double scale_factor = 0.95; 
    for (int i = 0; i < T1.size(); ++i) {
        double length = std::sqrt(T1[i] * T1[i] + T2[i] * T2[i]);

        if (length > 0) { 
            T1[i] /= length; 
            T2[i] /= length;
        }

        plt::plot({ 0.0, T1[i] * scale_factor }, { 0.0, T2[i] * scale_factor }, "orange");

        double arrow_size = 0.05;
        double angle = std::atan2(T2[i], T1[i]);

        double x1 = T1[i] - arrow_size * std::cos(angle - M_PI / 6);
        double y1 = T2[i] - arrow_size * std::sin(angle - M_PI / 6);
        double x2 = T1[i] - arrow_size * std::cos(angle + M_PI / 6);
        double y2 = T2[i] - arrow_size * std::sin(angle + M_PI / 6);

        plt::plot({ T1[i] * scale_factor, x1 * scale_factor }, { T2[i] * scale_factor, y1 * scale_factor }, "orange");
        plt::plot({ T1[i] * scale_factor, x2 * scale_factor }, { T2[i] * scale_factor, y2 * scale_factor }, "orange");

        plt::annotate(etiquetas[i], T1[i] * scale_factor, T2[i] * scale_factor);
    }

	//Grafico del circulo
    int num_points = 500; 
    std::vector<double> circle_x, circle_y;
    double radius = 1.0;  

    for (int i = 0; i < num_points; ++i) {
        double angle = 2 * 3.141592 * i / num_points;
        circle_x.push_back(radius * std::cos(angle));
        circle_y.push_back(radius * std::sin(angle));
    }

    //Datos y edicion del grafico
    plt::plot(circle_x, circle_y, "orange");  
    plt::xlabel(("Componente " + to_string(columna1) + " (" + to_string(x) + "%)"));
    plt::ylabel(("Componente " + to_string(columna2) + " (" + to_string(y) + "%)"));
    plt::title("Circulo de Correlacion");
    plt::axis("equal");
    plt::grid(true);
    plt::show();
}

int main(){

    cout << " == Matriz original ==" << endl;
    vector<string> nombres, materias;
    vector<vector<double>> original = leerCSV("EjemploEstudiantes.csv", nombres, materias);
    imprimirMatriz(original);
    cout << endl ;

    vector<vector<double>> estandarizada = estandarizar(original);
    cout << " == Matriz estandarizada ==" << endl;
    imprimirMatriz(estandarizada);
    cout << endl << endl;

    cout << " == Matriz de correlaciones/covarianzas ==" << endl;
    vector<vector<double>> correlacion = correlaciones(estandarizada);
    imprimirMatriz(correlacion);
    cout << endl << endl;

    pair<vector<double>, vector<vector<double>>> resultados = calcularValoresVectoresPropios(correlacion);

    vector<double> autovalores = resultados.first;
    vector<vector<double>> autovectores = resultados.second;

    cout << " == Valores propios ==" << endl;
    for (int i = 0; i < autovalores.size(); i++) {
        cout << autovalores[i] << " ";
    }
    cout << endl << endl;

    cout << " == Vectores propios ==" << endl;
    imprimirMatriz(autovectores);
    cout << endl ;

    cout << " == Matriz de componentes principales ==" << endl;
    vector<vector<double>> componentesP = componentes(estandarizada,autovectores);
    imprimirMatriz(componentesP);
    cout << endl ;

    cout << " == Matriz de calidades de individuos ==" << endl;
    vector<vector<double>> calidadesI = calidades(estandarizada,componentesP);
    imprimirMatriz(calidadesI);
    cout << endl;
  
    cout << " == Matriz de coordenadas de variables ==" << endl;
    vector<vector<double>> matrizT = coordenadasVariables(autovalores, autovectores);
    imprimirMatriz(matrizT);
    cout << endl;

    cout << " == Matriz de calidades de variables ==" << endl;
    vector<vector<double>> matrizS = calidadesVariables(matrizT, autovalores);
    imprimirMatriz(matrizS);
    cout << endl;

    cout << " == Vector de inercias ==" << endl;
    vector<double> vectorI = vectorInercias(autovalores);
    for (double i : vectorI) {
        cout << i << " ";
    }
    cout << endl;

    //Graficos
    PlanoPrincipal(nombres, componentesP, 0, 1);
    PlanoPrincipal(nombres, componentesP, 3, 4);
	Circulo_Correlacion(materias,matrizT,0,1, vectorI);
    Circulo_Correlacion(materias, matrizT, 3, 4, vectorI);

}

