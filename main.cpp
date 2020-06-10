#include <iostream>
#include <fstream>
#include <sstream>
#include <cassert>
#include <iomanip>
#include <cmath>
using namespace std;

class Matrix {
public:
    Matrix(); //konstruktor
    ~Matrix(); //destruktor
    void newMatrix(); //pamyat' pod matrizu
    void print() const; //pechat'
    int solve(); //opredelenie statusa resheniya: est' ili net
    int findNext(int *ar, int *array, int &place) const; //opredelenie nailuchshego varianta s naimen'shim kol-vom perestanovok
    double iteration(); //visheslenie iterazii i max delta
    int checkIsZeroes(); //proverka na 0
    int checkDUS(); //proverka DUS
    double getRowSumWithoutDiagonal(int row); //summa bez diagonal'nih elementov
    double getRowSum(int row); //summa elementov v stroke
    double ** getPermutatedRows(); //perestanovka
    double ** getPermutatedCoefficients(const int *permutations); //novaya matriza posle perestanovki
    int *getPermutations(int **swapsMatrix, int *swapsNumPerRow); //opredelenie nailuchshego varianta perestanovki (shema)
    double checkConvergence(); //proverka shodimosty
    void printSolutions(); //pechat' resheniya
    void init(const char* file); //chtenie iz faila i zapolnenie
private:
    double** matrix; //matriza
    double* solutions; //matriza resheniy
    int width; //shirina
    int height; //visota
    double epsilon; //epsilon
    static const int stepsNumForConvergence = 15; //shagi dlya iteraziy
    static const int stepsNumForDivergence = 5; //shagi dlya otslezhivaniya delta
};

Matrix::Matrix(){ //konstruktor
width = 0;
height = 0;
matrix = nullptr;
solutions = nullptr;
epsilon = 0;}

void Matrix::newMatrix() //pamyat' pod matrizu
{
    matrix = new double* [height];
    for (int i = 0; i < height; ++i)
        matrix[i] = new double[width];
}

void Matrix::init(const char* file) //chtenie iz faila i zapolnenie
{
    solutions = nullptr;
    ifstream in(file);
    in >> height >> width;
    in >> epsilon;
    newMatrix();
        for (int i = 0; i < height; ++i)
        for(int j = 0; j < width; ++j)
            in >> matrix[i][j];
    in.close();
}

int max_index(int *array, int size)	 //nomer max elementa
{
    int max = 0;
    for (int i = 0; i < size; i++) {
        if (array[i] > array[max])
            max = i;
    }
    return max;
}

Matrix::~Matrix()//destruktor
{
    if (matrix != nullptr)
        for (int i = 0; i < height; i++)
            delete[] matrix[i];
    delete[] solutions;
}

int Matrix::findNext(int *ar, int *array, int &place) const	//poisk varianta s naimen'shim kol-vom perestanovok
{
    int min = 0;
    place = -1;
    for (int i = 0; i < height; i++) {
        if (ar[i] == 0)
            return -1;
        if (ar[i] < ar[min]) {
            min = i;
            place = max_index(array + i * height, height);
        }
        else {
            if (ar[i] == ar[min]) {
                int max = max_index(array + i * height, height);
                if (array[i * height + max] > array[min * height + place]) {
                    min = i;
                    place = max;
                }
            }
        }
    }
    if (min == 0)
        place = max_index(array, height);
    return min;
}

int Matrix::checkIsZeroes(){ //proverka 0
    int isZeros = 0;
    for(int i = 0; i<height; i++)
        if (matrix[i][i]==0)
        {
            isZeros = 1;
            break;
        }
    return isZeros;
}

int Matrix::checkDUS(){ //proverka DUS
    int isAllDiagonalGTE = 0; //vse diagonalnie elementi bolshe ili ravni
    int isAnyDiagonalGT = 0; //hotyabi odin element strogo bolshe
    for(int i = 0; i < height; i++)
    {
        double sumWithoutDiagonal = getRowSumWithoutDiagonal(i);
        if (matrix[i][i] < sumWithoutDiagonal)
        {
            isAllDiagonalGTE = 0;
            break;
        }
        else
        {
            isAllDiagonalGTE = 1;
            if (matrix[i][i] > sumWithoutDiagonal)
            {
                isAnyDiagonalGT = 1;
            }
        }
    }
    return isAllDiagonalGTE && isAnyDiagonalGT;
}

double Matrix::getRowSumWithoutDiagonal(int row) //summa stroki bez diagonal'nogo elementa
{
    double sum = 0;
    for (int i = 0; i < height; ++i)
        if (i != row)
            sum += abs(matrix[row][i]);
    return sum;
}

double Matrix::getRowSum(int row) //summi strok
{
    double sum = 0;
    for (int i = 0; i < height; ++i)
        sum += abs(matrix[row][i]);
    return sum;
}

double ** Matrix::getPermutatedRows() //perestanovka
{
    int** swapsMatrix = new int*[height];
    for(int i=0; i < height; i++) {
        swapsMatrix[i] = new int[height];
    }
    int* swapsNumPerRow = new int[height];
    for(int i=0; i < height; i++) {
        swapsNumPerRow[i] = 0;
    }
    for(int i = 0; i<height; i++)
    {
        for(int j = 0; j<height; j++)
        {
            if (matrix[i][j] == 0)
            {
                swapsMatrix[j][i] = 0;
                continue;
            }
            swapsNumPerRow[i] ++;
            swapsMatrix[j][i] = 1;

            double otherSum = getRowSum(i) - abs(matrix[i][j]);
            if (abs(matrix[i][j]) == otherSum)
                swapsMatrix[j][i]++;
            if (abs(matrix[i][j]) > otherSum)
                swapsMatrix[j][i] += 2;
        }
    }

    int* permutations = getPermutations(swapsMatrix, swapsNumPerRow);

    for(int i = 0; i < height; ++i)
        delete[] swapsMatrix[i];
    delete[] swapsMatrix;
    delete[] swapsNumPerRow;

    if (permutations == nullptr)
        return nullptr;

    double** new_coefficients = getPermutatedCoefficients(permutations);
    delete[] permutations;

    return new_coefficients;
}

int *Matrix::getPermutations(int **swapsMatrix, int *swapsNumPerRow) { //shema perestanovki
    int* permutations = new int[height];
    for(int i = 0; i < height; ++i)
        permutations[i] = i;
    int where;
    int min;

    for (int k = 0; k < height; k++)
    {
        min = findNext(swapsNumPerRow, swapsMatrix[0], where); //poisk varianta s naimen'shim kol-vom perestanovok
        if (min == -1)
        {
            delete[] permutations;
            return nullptr;
        }
        for (int i = 0; i < height; i++)
            if (i != min && swapsMatrix[where][i])
            {
                swapsMatrix[where][i] = 0;
                swapsNumPerRow[i] -= 1;
            }
        permutations[min] = where;
        swapsNumPerRow[min] = 2*height + 1;
    }
    return permutations;
}

double ** Matrix::getPermutatedCoefficients(const int *permutations) { //novaya matriza posle perestanovki
    double** newCoefficients = new double*[height];
    for (int i = 0; i < height; i++)
        newCoefficients[i] = matrix[permutations[i]];
    return newCoefficients;
}

int Matrix::solve()	//opredelenie resheniya: est' ili net
{
    if (checkIsZeroes())
    {
        double ** permutatedCoefficients = getPermutatedRows();
        if (permutatedCoefficients == nullptr)
            return 2;
        else
        {
            delete[] matrix;
            matrix = permutatedCoefficients;
        }
    }
    int solvabilityStatus = !checkDUS();
    solutions = new double[height];
    for (int i = 0; i < height; i++)
        solutions[i] = 0;
    double delta = 0;
    if (solvabilityStatus==1)
    {
        delta = checkConvergence();
        if (delta < 0)
            return 1;
    }
    else
        delta = 2*epsilon;
    while (delta > epsilon) {
        delta = iteration();
    }
    return 0;
}

double Matrix::checkConvergence() //proverka shodimosty, otslezhivanie delta
{
    double prevDelta = 0;
    double delta = 0;
    int deltaDecreaseCounter = 0;
    for (int i = 0; i < stepsNumForConvergence; i++)
    {
        if (prevDelta > delta && abs(prevDelta - delta) >= epsilon)
            deltaDecreaseCounter++;
        else
            deltaDecreaseCounter = 0;
        prevDelta = delta;
        delta = iteration();
    }
    if (deltaDecreaseCounter < stepsNumForDivergence)
        delta = -1;
    return delta;
}

double Matrix::iteration()	//vicheslenie iterazii i max delta
{
    double maxDelta = 0;
    for (int i = 0; i < height; i++)
    {
        double prevSolutions = solutions[i];
        solutions[i] = matrix[i][height];
        for (int j = 0; j < height; j++)
            if (j != i)
                solutions[i] -= solutions[j] * matrix[i][j];

        solutions[i] /= matrix[i][i];
        if (abs(solutions[i]) < epsilon)
            solutions[i] = 0;
        prevSolutions = abs(prevSolutions - solutions[i]);
        maxDelta = max(prevSolutions, maxDelta);
    }
    return maxDelta;
}

void Matrix::print() const	//prints the system (matrix only)
{
    cout << "MATRIX: " << endl;
    for (int i = 0; i < height; i++) {
        for (int j = 0; j < width - 1; j++)
            cout << setw(15) << scientific << matrix[i][j] << " ";
        cout << setw(15) << scientific << matrix[i][width - 1] << " " << endl;
    }
    cout << endl;
}

void Matrix::printSolutions() //pechat' resheniy
{
    init("system.txt");
    print();
    int solvabilityStatus;
    if (solutions == nullptr)
         solvabilityStatus = solve();
    switch (solvabilityStatus)
    {
        case 0:
        {
            cout <<  "SOLUTIONS:" << endl;
            for (int i = 0; i < height; i++)
                cout  << solutions[i] << " " << endl;
            break;
        }
        case 1:
        {
            print();
            cout << "no solutions: DUS isn't executed & matrix diverges";
            break;
        }
        case 2:
        {
            print();
            cout << "no solutions: zeroes exist & no appropriate change";
        }
        default:
        {
            cout << "unknown error";
            break;
        }

    }

}

int main()
{
    Matrix system;
    system.printSolutions();
}
