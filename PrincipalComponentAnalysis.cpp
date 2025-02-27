#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include <string>
#include <cmath>
#include <iterator>
#include <Eigen/Eigenvalues>
#include "matplotlibcpp.h"

using namespace std;
using namespace Eigen;
namespace plt = matplotlibcpp;

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

vector<vector<double>> leerCSV(const string& nombreArchivo) {
    ifstream archivo(nombreArchivo);
    vector<vector<double>> datos;
    string linea;

    if (!archivo.is_open()) {
        cerr << "Error al abrir el archivo" << endl;
        return datos;
    }

    getline(archivo, linea);
    while (getline(archivo, linea)) {
        stringstream ss(linea);
        string valor;
        vector<double> fila;

        // skipea los nombres
        getline(ss, valor, ';');

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
		cout << avg << endl;

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

pair<vector<double>, vector<vector<double>>> calcularValoresVectoresPropios(const vector<vector<double>>& matriz) {
    int row = matriz.size();  

    // convertir de tipo vector a MatrixXd
    MatrixXd mat(row, row);
    for (int i = 0; i < row; i++) {
        for (int j = 0; j < row; j++) {
            mat(i, j) = matriz[i][j];
        }
    }

    SelfAdjointEigenSolver<MatrixXd> solver(mat);

    VectorXd valoresPropios = solver.eigenvalues();  
    MatrixXd vectoresPropios = solver.eigenvectors();


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

vector<vector<double>> calidades(const vector<vector<double>>& matrizX, const vector<vector<double>>& matrizC) {
    vector<vector<double>> matrizQ;
    vector<double> temp;

    int row = matrizX.size();
    int col = matrizX[0].size();
    double calculo = 0, suma = 0;//variable para los calculos de cada elemento de la matriz

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

int main(){

    cout << " == Matriz original ==" << endl;
	vector<vector<double>> original = leerCSV("EjemploEstudiantes.csv");
    imprimirMatriz(original);
    cout << endl << endl;

    vector<vector<double>> estandarizada = estandarizar(leerCSV("EjemploEstudiantes.csv"));
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

    plt::plot({ 1, 2, 3, 4 }, "*");
    plt::show();
    plt::detail::_interpreter::kill();

}

