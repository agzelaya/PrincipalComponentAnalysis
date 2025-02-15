#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include <string>
#include <cmath>

using namespace std;

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


void imprimirMatriz(const vector<vector<double>>& matriz) {
    for (const auto& fila : matriz) {
        for (double valor : fila) {
            cout << valor << "\t";
        }
        cout << endl;
    }
}

int main(){
    imprimirMatriz(leerCSV("EjemploEstudiantes.csv"));
}

