
#include <C_General.hpp>
#include <C_File.hpp>
#include <C_Arguments.hpp>
#include <C_Trace.hpp>
#include <C_Image.hpp>
#include <iostream>
#include <stdio.h>
#include <string>
#include <algorithm>
#include <iterator>
#include <map>
#include <vector>
#include <cmath>
#include <conio.h>
#include <limits>
#include <thread>         //sleep
#include <chrono>  

//Para la configuracion de la barra de progreso
#define PBSTR "||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||"
#define PBWIDTH 60

//DECLARACION CAMPOS GLOBALES
C_Image image1, image2, image3;
C_Matrix::IndexT firstRow1, lastRow1, firstCol1, lastCol1, rowN1, colN1;
C_Matrix::IndexT firstRow2, lastRow2, firstCol2, lastCol2, rowN2, colN2;
C_Matrix::IndexT firstRow3, lastRow3, firstCol3, lastCol3, rowN3, colN3;
string message;

//Si hay una correlacion perfecta (100%), no hace falta suavizar la linea
bool exacto = false;


// Mapas donde se guardaran los valores de Corr.Cruzada (KEYS) y las posiciones espaciales de las subimagenes 1 y 2
// de la imagen original (VALUES)
map <double, vector<long>> valoresCorrelacion1;
map <double, vector<long>> valoresCorrelacion2;


//INICIALIZACION DE CAMPOS Y FUNCIONES
void Inicializacion_Variables() {

	//Primeras y ultimas filas de cada imagen
	firstRow1 = image1.FirstRow();
	lastRow1 = image1.LastRow();
	firstRow2 = image2.FirstRow();
	lastRow2 = image2.LastRow();
	firstRow3 = image3.FirstRow();
	lastRow3 = image3.LastRow();

	//Primeras y ultimas columnas de cada imagen
	firstCol1 = image1.FirstCol();
	lastCol1 = image1.LastCol();
	firstCol2 = image2.FirstCol();
	lastCol2 = image2.LastCol();
	firstCol3 = image3.FirstCol();
	lastCol3 = image3.LastCol();

	//Numero de filas y columnas de cada imagen
	rowN1 = image1.RowN();
	colN1 = image1.ColN();
	rowN2 = image2.RowN();
	colN2 = image2.ColN();
	rowN3 = image3.RowN();
	colN3 = image3.ColN();
	
	//String para mostrar por consola
	message = "";

}


//Para mostrar por pantalla el progreso del algoritmo. 
string Mini_Consola(string mensaje) {

	cout << mensaje;

	return mensaje;
}



//Limite de anchura
long PorcentajeAnchura(long long numero) {

	long long porcentaje = (numero*10)/100;

	return porcentaje;
}


//Limite de altura
long PorcentajeAltura(long long numero) {

	long long porcentaje = (numero * 10) / 100;

	return porcentaje;
}


//Se realiza el calculo del indice de correlacion cruzada, entre dos submatrices de las imagenes originales
double Correlacion_Cruzada(C_Image mat1, C_Image mat2, int k) {

	C_Matrix aux1 = C_Matrix(mat1);
	C_Matrix aux2 = C_Matrix(mat2);
	C_Matrix aux3 = C_Matrix(mat2);
	double correlation;
	long long int sum1, sum2, sum3;
	
	aux1.SubtractEscalar(aux1.Mean());
	aux2.SubtractEscalar(aux2.Mean());

	aux3.MultiplyElm(aux1, aux2);

	sum1 = aux3.Sum();

	aux3.MultiplyElm(aux1, aux1);
	sum2 = aux3.Sum();

	aux3.MultiplyElm(aux2, aux2);
	sum3 = aux3.Sum();
	
	long long int temp1 = sum2 * sum3;
	correlation = sum1 / sqrt(sum2 * sum3);

	return correlation;
}


//Mostrar por pantalla el porcentaje del proceso
void barraDeProgreso(double porcentaje) {
	int val = (int)(porcentaje * 100);
	int lpad = (int)(porcentaje * PBWIDTH);
	int rpad = PBWIDTH - lpad;
	printf("\r%3d%% [%.*s%*s]", val, lpad, PBSTR, rpad, "");
	fflush(stdout);
}



//COMPARACION DE MATRICES (IMAGENES) MEDIANTE LA CORRELACION CRUZADA 
void Comparacion_De_Imagenes(C_Matrix mat1, C_Matrix mat2) {

	valoresCorrelacion1.clear();
	valoresCorrelacion2.clear();
	long SubColBeginning1 = mat1.ColN() / 2;									//Subimagen1 (rectangulo central): columna inicial y   
	long SubColEnd = ((mat1.ColN() / 2) + PorcentajeAnchura(mat1.ColN() / 2));	//columna final
	long SubColBeginning2 = mat1.FirstCol();
	long SubColEnd2 = (SubColBeginning2 + PorcentajeAnchura(mat1.ColN() / 2));
	long RowBeginning = mat1.FirstRow();	//Asumiendo que ambas imagenes poseen la misma resolucion y por lo tanto, 						
	long RowEnd = mat1.LastRow();			//mismo numero de columnas
	long actualRowSub2 = RowBeginning;
	double correlacion, correlacionAux;
	bool max = false;			//El porcentaje es 100% 
	bool principio = false;			
	double porcentaje;
	int m = 0;

	//Submatriz 1: Extremo izquierdo de la imagen 1
	C_Image subMat1 = C_Matrix(mat1, mat1.FirstRow(), mat1.LastRow(), SubColBeginning2, SubColEnd2, mat1.FirstRow(), SubColBeginning2);
	
	//Submatriz 1: Mitad de la imagen 1
	C_Image temp1 = C_Matrix(mat1, mat1.FirstRow(), mat1.LastRow(), SubColBeginning1, SubColEnd, mat1.FirstRow(), SubColBeginning1);
	
	//Submatriz 2
	C_Image subMat2 = C_Matrix(mat2, mat2.FirstRow(), subMat1.LastRow(), mat2.FirstCol(), subMat1.ColN(), mat2.FirstRow(), mat2.FirstCol());

	correlacionAux = Correlacion_Cruzada(subMat1, subMat2, m);
	correlacion = Correlacion_Cruzada(temp1, subMat2, m);
	long limit1 = mat1.ColN() - subMat1.ColN();
	long limit2 = (mat1.ColN() / 2) - subMat1.ColN();
	long limit = limit1;

	//Si empezamos desde la mitad de la img1 o desde el principio de la img1
	if (correlacionAux > correlacion) {
		principio = true;
		temp1 = subMat1;
		SubColEnd = SubColEnd2;
		limit = limit2;
	}
	
	//Mover subimagen1 de izquierda a derecha
	for (int i = SubColEnd; i < limit; i++) { //ultima columna de la subimg1, evitar que salga de la img1

		porcentaje = (double)i / (double)limit;
		barraDeProgreso(porcentaje);

		if (max) break;

		//Mover subimagen1 hacia "arriba" (recortamos imagen)
		for (int j = RowEnd; j > (RowEnd - PorcentajeAltura(RowEnd)); j--) {

			if (max) break;

			C_Image temp1 = C_Matrix(mat1, mat1.FirstRow(), j, 
				(principio ? SubColBeginning2 : SubColBeginning1), i, mat1.FirstRow(), 
				(principio ? SubColBeginning2 : SubColBeginning1));
			
			long actualRowSub2 = RowBeginning;

			//Ahora se mueve de arriba a abajo la subimg2 (considerar todas las posibilidades)
			for (int k = j; k <= RowEnd; k++) {

				//La subimg2 debe tener las mismas dimensiones que la subimg1
				C_Image temp2 = C_Matrix(mat2, actualRowSub2, k, mat2.FirstCol(), temp1.ColN(), 
										actualRowSub2, PorcentajeAnchura(colN2));

				correlacion = Correlacion_Cruzada(temp1, temp2, k);

					//Si hemos encontrado el punto, paramos de buscar.
					if (correlacion == 1) {
						max = true;
						exacto = true;
						vector<long> arr1;
						arr1.push_back(temp1.FirstRow());
						arr1.push_back(temp1.LastRow());
						arr1.push_back(temp1.FirstCol());
						arr1.push_back(temp1.LastCol());
						vector<long> arr2;
						arr2.push_back(temp2.FirstRow());
						arr2.push_back(temp2.LastRow());
						arr2.push_back(temp2.FirstCol());
						arr2.push_back(temp2.LastCol());

						valoresCorrelacion1.emplace(correlacion, arr1);
						valoresCorrelacion2.emplace(correlacion, arr2);
						break;
					}

				//0 -> firstRow, 1 -> lastRow, 2 -> firstCol, 3 -> lastCol
				vector<long> arr1;
				arr1.push_back(temp1.FirstRow());
				arr1.push_back(temp1.LastRow());
				arr1.push_back(temp1.FirstCol());
				arr1.push_back(temp1.LastCol());
				vector<long> arr2;
				arr2.push_back(temp2.FirstRow());
				arr2.push_back(temp2.LastRow());
				arr2.push_back(temp2.FirstCol());
				arr2.push_back(temp2.LastCol());

				valoresCorrelacion1.emplace(correlacion, arr1);
				valoresCorrelacion2.emplace(correlacion, arr2);

				actualRowSub2++;
			}
		}

		//Mover subimagen1 hacia "abajo"
		for (int j = RowBeginning; j < PorcentajeAltura(subMat2.RowN()); j++) {

			if (max) break;

			C_Image temp1 = C_Matrix(mat1, j, RowEnd, (principio ? SubColBeginning2 : SubColBeginning1), i, mat1.FirstRow(), 
													  (principio ? SubColBeginning2 : SubColBeginning1));
			long actualRowSub2 = RowBeginning; //subimg2 siempre desde el comienzo de la imagen original

			for (int k = (RowEnd + 1) - j; k <= RowEnd; k++) {

				//Ahora se mueve de arriba a abajo la subimg2 (considerar todas las posibilidades)
				C_Image temp2 = C_Matrix(mat2, actualRowSub2, k, mat2.FirstCol(), temp1.ColN(), 
										 actualRowSub2, PorcentajeAnchura(mat2.ColN()));

				correlacion = Correlacion_Cruzada(temp1, temp2, k);

					//Si hemos encontrado el punto, paramos de buscar.
					if (correlacion == 1) {
						max = true;
						exacto = true;
						vector<long> arr1;
						arr1.push_back(temp1.FirstRow());
						arr1.push_back(temp1.LastRow());
						arr1.push_back(temp1.FirstCol());
						arr1.push_back(temp1.LastCol());
						vector<long> arr2;
						arr2.push_back(temp2.FirstRow());
						arr2.push_back(temp2.LastRow());
						arr2.push_back(temp2.FirstCol());
						arr2.push_back(temp2.LastCol());

						valoresCorrelacion1.emplace(correlacion, arr1);
						valoresCorrelacion2.emplace(correlacion, arr2);
						break;
					}
				

				//0 -> firstRow, 1 -> lastRow, 2 -> firstCol, 3 -> lastCol
				vector<long> arr1;
				arr1.push_back(temp1.FirstRow());
				arr1.push_back(temp1.LastRow());
				arr1.push_back(temp1.FirstCol());
				arr1.push_back(temp1.LastCol());
				vector<long> arr2;
				arr2.push_back(temp2.FirstRow());
				arr2.push_back(temp2.LastRow());
				arr2.push_back(temp2.FirstCol());
				arr2.push_back(temp2.LastCol());

				valoresCorrelacion1.emplace(correlacion, arr1);
				valoresCorrelacion2.emplace(correlacion, arr2);

				actualRowSub2++;
			}
		}
		if (!max) (principio ? SubColBeginning2++ : SubColBeginning1++);
	}
}


//Suavizar la linea de separacion de las dos imagenes para que se fusionen correctamente
C_Image SuavizadoLinea(C_Matrix imageFinal, C_Matrix image1, C_Matrix image2, 
						vector<long> directions1, vector<long> directions2, long puntoMedio) {


	long anchura = ((long)((directions1.at(3) - directions1.at(2))*3)/4); //Anchura de la "submatriz" en la cual se aplicara el suavizado
	anchura = (long)((abs(anchura))/2);									  //La anchura debe ser un valor positivo
	long inicio = abs((anchura/2) - puntoMedio) + 1;
	long fin = abs((anchura / 2) + puntoMedio) + 1;
	double porcentajeIMG1 = 1;
	double porcentajeIMG2 = 0;
	
	long firstColsub1 = directions1.at(2);
	long lastColsub1 = directions1.at(3);
	
	long firstRowsub1 = directions1.at(0);
	long firstRowsub1AUX = firstRowsub1;
	long lastRowsub1 = directions1.at(1);
	
	long firstColsub2 = directions2.at(2);
	long lastColsub2 = directions2.at(3);
	
	long firstRowsub2 = directions2.at(0);
	long firstRowsub2AUX = firstRowsub2;
	long lastRowsub2 = directions2.at(1);
	
	double porcentaje = 1/(double)(fin - inicio);


	for (int j = inicio; j <= fin; j++) {
		firstRowsub1 = firstRowsub1AUX;
		firstRowsub2 = firstRowsub2AUX;

		for (int i = imageFinal.FirstRow(); i <= imageFinal.LastRow(); i++) {
			
			imageFinal(i, j) =	abs((porcentajeIMG1 * image1(firstRowsub1, firstColsub1)) +
								   ((porcentajeIMG2 * image2(firstRowsub2, firstColsub2))));
			
			if (firstRowsub1 < lastRowsub1)
				firstRowsub1++;
			if(firstRowsub2 < lastRowsub2)
				firstRowsub2++;
		}
		
		if(firstColsub1 < lastColsub1)
			firstColsub1++;
		if (firstColsub2 < lastColsub2)
			firstColsub2++;
		
		porcentajeIMG1 -= porcentaje;
		porcentajeIMG2 += porcentaje;
	}

	return imageFinal;
}




//UNIR IMAGENES (IMG1 con IMG2, y luego IMG12 con IMG3)
C_Image Unir_Imagenes(C_Image img1, C_Image img2) {

	//0 -> firstRow, 1 -> lastRow, 2 -> firstCol, 3 -> lastCol
	vector<long> finalDirections1 = (valoresCorrelacion1.rbegin())->second; //Mapas se ordenan de forma ascendente por sus llaves (keys)
	vector<long> finalDirections2 = (valoresCorrelacion2.rbegin())->second; //Por lo tanto los valores de la última key
																			//(mayor valor de correlacion) son los necesarios

	//Determinamos la altura que tendra la imagen combinada
	long lastRowFinal = (finalDirections1.at(1) < finalDirections2.at(1)) ? finalDirections1.at(1) : finalDirections2.at(1);
	long firstRowFinal = (finalDirections1.at(0) > finalDirections2.at(0)) ? finalDirections1.at(0) : finalDirections2.at(0);

	long dimensionsub2 = (finalDirections2.at(3) - finalDirections2.at(2)); 
	dimensionsub2 += dimensionsub2;										   //IMPORTANTE: Mover el indice para que coincida correctamente 
	long lastColFinal = img1.LastCol() + (img2.LastCol() - dimensionsub2); //Es posible que se tenga que recortar más adelante
	long limiteimg1 = finalDirections1.at(2) - dimensionsub2;			   //IMPORTANTE DONDE SE CORTA LA PRIMERA IMAGEN

	//Imagen "plantilla" que posee las dimensiones finales y se rellenara con el contenido de las imagenes 1 y 2
	C_Image combinada12 = C_Image(1, lastRowFinal, img1.FirstCol(), lastColFinal, 0);

	int k = 1;
	//Rellenamos la imagen plantilla con la img1
	for (int i = finalDirections1.at(0); i <= finalDirections1.at(1); i++) {

		
		for (int j = img1.FirstCol(); j <= limiteimg1; j++) {

			combinada12(k, j) = img1(i, j);		
		}

		if(k < combinada12.RowN())
		k++;
	}
		
	k = 1;
	//Luego rellenamos combinada12 con la img2 (finalDirections1 firstCol + img2.firstCol())
	for (int i = finalDirections2.at(0); i <= finalDirections2.at(1); i++) {
		
		for (int j = img2.FirstCol(); j <= img2.LastCol(); j++) {

			//Menos 1 para tener en cuenta la primera columna tambien
			combinada12(k, (limiteimg1 + j)) = img2(i, j);
		}
		if (k < combinada12.RowN())
		k++;
	}
	 
	//Recortamos el espacio sobrante (eliminacion de franjas negras horizontales)
	long limiteFinal = limiteimg1 + img2.LastCol();
	
	//IMAGEN FINAL
	C_Image imgFinal = C_Matrix(combinada12, combinada12.FirstRow(), combinada12.LastRow(), 
		                        combinada12.FirstCol(), limiteFinal, combinada12.FirstRow(), combinada12.FirstCol());

	//Si no se ha encontrado el punto exacto donde juntar las imagenes, se suaviza la linea de separacion de ambas
	if(exacto == false)
	imgFinal = SuavizadoLinea(imgFinal, img1, img2, finalDirections1, finalDirections2, limiteimg1);

	return imgFinal;
}


int  main(int argc, char **argv) {

	C_Print("ALGORITMO: IMAGEN PANORAMICA");
	C_Print("-----------------------------");

	char nombreImagen1[50];
	char nombreImagen2[50];
	char nombreImagen3[50];
	char numImg;


	C_Print("Cuantas imagenes va a unir? (Elija 2 o 3)");
	cin >> numImg;
	cin.ignore();

	C_Print("Por favor, introduzca el nombre de la imagen 1 y su extension");
	cin.getline(nombreImagen1, 50);

	C_Print("Por favor, introduzca el nombre de la imagen 2 y su extension");
	cin.getline(nombreImagen2, 50);

	if(numImg == '3') {
	C_Print("Por favor, introduzca el nombre de la imagen 3 y su extension");
	cin.getline(nombreImagen3, 50);
	}

	//Lectura de imagenes
	message = Mini_Consola("Leyendo imagenes...");
	std::this_thread::sleep_for(std::chrono::seconds(2));
	image1.ReadBMP(nombreImagen1);
	image2.ReadBMP(nombreImagen2);
	if (numImg == '3')
	image3.ReadBMP(nombreImagen3);

	int error = 0;
	//Diferentes Fail() con el objetivo de que si alguna de las imagenes tiene un problema, podemos saber de cual se trata.
	C_IfError(image1.Fail(), "Error: No existe la imagen 1, ha introducido mal el nombre o formato incorrecto.", error++);
	C_IfError(image2.Fail(), "Error: No existe la imagen 2, ha introducido mal el nombre o formato incorrecto.", error++);
	C_IfError(image3.Fail(), "Error: No existe la imagen 3, ha introducido mal el nombre o formato incorrecto.", error++);

	if (error > 0) {
		C_Print("Reinicie el programa e intentelo de nuevo.");
		system("pause");
		return -1;
	}
	message += Mini_Consola("OK\n");	//Lectura correcta
	

	//Inicializacion de variables
	message = Mini_Consola("Inicializando variables...");
	std::this_thread::sleep_for(std::chrono::seconds(2));
	Inicializacion_Variables();
	message += Mini_Consola("OK\n");	
	

    //Conversion de imagenes en escala de grises
	message = Mini_Consola("Convirtiendo imagenes en escala de grises...");
	std::this_thread::sleep_for(std::chrono::seconds(2));
	image1.Grey();
	image2.Grey();
	if (numImg == '3')
	image3.Grey();
	message += Mini_Consola("OK\n");
	
	//Desarrollo de la comparacion y union para 2 imagenes
	C_Print(" ");
	C_Print("**************************");
	C_Print("Analizando imagenes 1 y 2");
	C_Print("**************************");
	C_Print("Dependiendo del tamanio de la imagen este proceso podria tardar unos minutos, espere por favor...");
	Comparacion_De_Imagenes(image1, image2);
	C_Print(" ");
	message = Mini_Consola("Uniendo imagenes...");
	C_Image unida1 = Unir_Imagenes(image1, image2);
	message += Mini_Consola("OK\n");
	C_Print(" ");

	if(exacto == false)
		C_Print("Suavizando linea de separacion entre imagenes...");

	//Desarrollo de la comparacion y union para 3 imagenes
	if (numImg == '3') {

		C_Print(" ");
		C_Print("**************************");
		C_Print("Analizando imagenes 2 y 3");
		C_Print("**************************");
		C_Print(" ")
		Comparacion_De_Imagenes(unida1, image3);
		C_Print(" ");
		message = Mini_Consola("Uniendo imagenes...");
		C_Image imageFinal = Unir_Imagenes(unida1, image3);
		message += Mini_Consola("OK\n");
		
		if (exacto == false)
			C_Print("Suavizando linea de separacion entre imagenes...");
		
		imageFinal.WriteBMP("Imagen_Unida.bmp");
		message += Mini_Consola("OK\n");
		C_Print(" ");
	} else
		unida1.WriteBMP("Imagen_Unida.bmp");

	C_Print("Busca la imagen: Imagen_Unida.bmp en la carpeta Run ");
	C_Print(" ");
	C_Print("**************************************");
	C_Print("IMAGEN UNIDA");
	C_Print("Programa Finalizado");
	C_Print("Autor: Angel Ciudad Montalban");
	C_Print("Copyright (c) 2019");
	C_Print("**************************************");
	C_Print(" ")
	system("pause");
}

//return Test(argc, argv);

int Test(int argc, char **argv);
