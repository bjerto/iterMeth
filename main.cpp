#include <iostream>
#include <fstream>
#include <sstream>
#include <cassert>
#include <iomanip>
#include <cmath>
using namespace std;

class Matrix {
public:
    Matrix();
    Matrix(const char* file);
    ~Matrix();
    void print() const; //prints the system (coefficients only)
    void prn_solutions(bool sci = true); //prints the solutions only
private:
    void sub(int from, int row, double val); //subtracts the row "row" from the row "from"
    void swap(double*& a, double*& b); //swap the rows coefficients and b
    void error(int err) const; // prints an error based on coefficients code
    void rearrange(int n[]);	//rearranges the Matrix according to coefficients given permutation
    double* sums() const; //calculates the sums of elements
    int solve(); //calculate the answers; errors: 0 - the system is solved; 1 - the system is solvable but the row diverges; 2 - the system is not get_solvability_status because of zeros
    int getSolvabilityStatus(int &isZeros, int *permutation); //return 0 - is 100% solvable, 1 - solve under control, 2 - totally unsolvable, replaces if 0's are met;
    int nonzero(int n) const; // finds the first string with non-zero element
    int find_next(int n[], int b[], int& place) const; //returns id of the next element to choose place for and where you should place it
    int find_max(int* b, int sz) const; //returns id of the largest element in array coefficients.
    int zero(double n) const {	return (abs(n) < epsilon);	}	//checks if n tends to zero
    double step(); //calculates one step and returns maximum delta
    double** coefficients;
    double* solutions; 
    int width;
    int height;
    int checkIsZeros();
    int checkDUS();
    double rowSumWithoutDiagonal(int row);
    double epsilon; //used for checking if the row diverges
    static const int delay = 15, wait = 5; //amount of steps done without checking for divergence and amount of steps allowing delta to increase
};

Matrix::Matrix(): width(0), height(0), coefficients(nullptr), solutions(nullptr), epsilon(0) {}

Matrix::Matrix(const char* file): solutions(nullptr)
{
    ifstream in(file);
    in >> height >> width;
    in >> epsilon;

    coefficients = new double* [height];
    for (int i = 0; i < height; ++i)
        coefficients[i] = new double[width];

    for (int i = 0; i < height; ++i)
        for(int j = 0; j < width; ++j)
            in >> coefficients[i][j];
    in.close();
}

Matrix::~Matrix()
{
    if (coefficients != nullptr)
        for (int i = 0; i < height; i++)
            delete[] coefficients[i];
    delete[] solutions;
}






int Matrix::find_next(int n[], int b[], int& place) const	//finds coefficients string with the least amount of solutions and the best solution
{
    int min = 0; //id of coefficients row with min amount of possible replacements
    place = -1; //coefficients place to put our string to
    for (int i = 0; i < height; i++) {
        if (n[i] == 0)
            return -1; //if some row has no place to put it, this system is not solvable
        if (n[i] < n[min]) {
            min = i; //if some row has less solutions than the current one then
            place = find_max(b + i * height, height); //recalculate the best solution place
        }
        else {
            if (n[i] == n[min]) { //if they have the same amount of solutions
                int max = find_max(b + i * height, height); //finding the best solution for coefficients row i
                if (b[i * height + max] > b[min * height + place]) { //if the best solution for i is better than the best solution for current min then replace
                    min = i;
                    place = max;
                }
            }
        }
    }
    if (min == 0)
        place = find_max(b, height);
    return min;
}

int Matrix::find_max(int* b, int sz) const	 //returns id of the largest element in array coefficients.
{
    int max = 0;
    for (int i = 0; i < sz; i++) {
        if (b[i] > b[max])
            max = i;
    }
    return max;
}

void Matrix::sub(int from, int row, double val)	//subtracts the row "row" from the row "from"
{
    for (int i = row; i < width; i++)
        coefficients[from][i] -= coefficients[row][i] * val;
}

void Matrix::swap(double*& a, double*& b)	//swap the rows coefficients and b
{
    double* temp = a;
    a = b;
    b = temp;
}

void Matrix::rearrange(int n[])	//rearranges the Matrix according to coefficients given permutation
{
    double** temp = new double*[height];
    for (int i = 0; i < height; i++)
        temp[i] = coefficients[n[i]];
    delete[] coefficients;
    coefficients = temp;
}

int Matrix::nonzero(int n) const	//returns the first non-zero element met
{
    int i;
    for (i = n; i < height; i++)
        if (!zero(coefficients[i][n]))
            break;
    return i;
}

double* Matrix::sums() const	//calculates the sums of elements
{
    double* sum = new double[height];
    for (int i = 0; i < height; i++) {
        sum[i] = 0;
        for (int j = 0; j < height; j++)
            sum[i] += abs(coefficients[i][j]);
    }
    return sum;
}

int Matrix::checkIsZeros(){
    int isZeros = 0; //shows if there are zero diagonal elements in the original Matrix => if you should replace things
    for(int i = 0; i<height; i++)	//making coefficients swap table. 0 - swap is impossible (0 on diagonal), 1 - swap is possible, but the requirement is not met, 2 - min requirement is met, 3 - max requirement is met;
        if (coefficients[i][i]==0)
        {
            isZeros = 1;	//if there is coefficients 0 on coefficients diagonal element, raise zero flag and reset the rest (we will be swapping it)
            break;
        }
    return isZeros;
}

int Matrix::checkDUS(){
    int isAllDiagonalGTE = 0; //shows if all of the diagonal elements are at least equal to the sum of others
    int isAnyDiagonalGT = 0; //shows if at least one diagonal element in any row is strictly larger than the sum of absolute amounts of others
    double* sum = sums(); //calculating sums of elements and store them
    for(int i = 0; i < height; i++)	//making coefficients swap table. 0 - swap is impossible (0 on diagonal), 1 - swap is possible, but the requirement is not met, 2 - min requirement is met, 3 - max requirement is met;
    {
        double sumWithoutDiagonal = rowSumWithoutDiagonal(i);	//checking which swaps would do
        if (coefficients[i][i] < sumWithoutDiagonal)
        {
            isAllDiagonalGTE = 0;
            break;
        }
        else
        {
            isAllDiagonalGTE = 1;
            if (coefficients[i][i] > sumWithoutDiagonal)
            {
                isAnyDiagonalGT = 1;
            }
        }
    }
    delete[] sum;
    return isAllDiagonalGTE && isAnyDiagonalGT;
}

double Matrix::rowSumWithoutDiagonal(int row)
{
    double sum = 0;
    for (int i = 0; i < height; ++i)
        if (i != row)
            sum += abs(coefficients[row][i]);
    return sum;
}

int Matrix::getSolvabilityStatus(int &isZeros, int *permutation){ //return 0 - is 100% solution, 1 - solve under control, 2 - totally unsolvable, replaces if 0's are met;
    isZeros = checkIsZeros();//shows if there are zero diagonal elements in the original Matrix => if you should replace things
    if (!isZeros)
    {
        return !checkDUS();
    }
    else
    {
        permutateRows();
        int isAllDiagonalGTE = 1; //shows if all of the diagonal elements are at least equal to the sum of others
        double* sum = sums(); //calculating sums of elements and store them
        int swapsMatrix[height][height]; //Matrix of swaps that can be done to make the system solvability status
        int* swapsNumPerRow = new int[height]; //number of possible swaps for the row
        for(int i=0; i < height; i++) {
            swapsNumPerRow[i] = 0;
        }
        int isAnyDiagonalGT = 0; //shows if at least one diagonal element in any row is strictly larger than the sum of absolute amounts of others
        for(int i = 0; i<height; i++)	//making coefficients swap table. 0 - swap is impossible (0 on diagonal), 1 - swap is possible, but the requirement is not met, 2 - min requirement is met, 3 - max requirement is met;
        {
            for(int j = 0; j<height; j++){
                double oth = sum[i]-abs(coefficients[i][j]);	//checking which swaps would do
                if (coefficients[i][j]!=0){
                    swapsNumPerRow[i] ++;
                    swapsMatrix[j][i] = 1;
                    if (abs(coefficients[i][j]) >= oth){
                        swapsMatrix[j][i]++;
                        if (abs(coefficients[i][j]) > oth){
                            swapsMatrix[j][i]++;
                            if (i==j) {
                                isAnyDiagonalGT = 1;
                            }
                        }
                    }
                    else {
                        if (i==j) {
                            isAllDiagonalGTE = 0;
                        }
                    }
                }
                else{
                    swapsMatrix[j][i] = 0;
                    if (i==j){
                        isAllDiagonalGTE = 1;
                        isAnyDiagonalGT = 0;
                    }
                }
            }
        }
        delete[] sum;

        int where;
        int min;

        for (int k = 0; k < height; k++) {
            min = find_next(swapsNumPerRow, swapsMatrix[0], where); //find coefficients row with the least amount of options
            if (min == -1) //if there are no non-zero options, leave with error code
                return 2;
            if (swapsMatrix[min][where] == 1) //if it is not even equal to the sum, zero the flag, we'll have to trace
                isAllDiagonalGTE = 0;
            else if (swapsMatrix[min][where] == 3) //if the diagonal element is larger than the sum, raise the flag
                isAnyDiagonalGT = 1;
            for (int i = 0; i < height; i++)	//cross it out for other rows
                if (i != min)
                    if (swapsMatrix[where][i]) {
                        swapsMatrix[where][i] = 0;
                        swapsNumPerRow[i] -= 1;
                    }
            permutation[min] = where;	//save this permutation
            swapsNumPerRow[min] = 2*height + 1;	//so that we wont come back to this row
        }
    return !(isAllDiagonalGTE && isAnyDiagonalGT);
    }
}

int Matrix::solve()	//calculate the answers; errors: 0 - the system is solved; 1 - the system is solvable but the row diverges; 2 - the system is not get_solvability_status because of zeros
{
    int permutations[height];
    int isZeros = 0;
    int solvabilityStatus = getSolvabilityStatus(isZeros, permutations);
    if (solvabilityStatus == 2)
        return 2;
    if (isZeros)
        rearrange(permutations);	//use our new permutation
    print();
    solutions = new double[height];
    for (int i = 0; i < height; i++)
        solutions[i] = 0;
    double prev_delta = 0;
    double delta = 0;
    int more = 0;
    if (solvabilityStatus==1){
        cout<<"Making "<<delay<<" trial steps to make sure the row does not diverge."<<endl;
        for (int i = 0; i < delay; i++) {
            if (prev_delta > delta && !zero(prev_delta - delta))
                more++;
            else
                more = 0;
            prev_delta = delta;
            delta = step();
        }
        if (more < wait)
            return 1;
    }
    else
        delta = 2*epsilon;
    while (delta > epsilon) {
        delta = step();
    }
    return 0;
}

double Matrix::step()	//calculates one step and returns maximum delta
{
    double max_delta = 0;
    for (int i = 0; i < height; i++) {
        double temp = solutions[i];
        solutions[i] = coefficients[i][height];
        for (int j = 0; j < height; j++) {
            if (j == i)
                continue;
            solutions[i] -= solutions[j] * coefficients[i][j];
        }
        solutions[i] /= coefficients[i][i];
        if (zero(solutions[i]))
            solutions[i] = 0;
        temp = abs(temp - solutions[i]);
        max_delta = max(temp, max_delta);
    }
    return max_delta;
}

void Matrix::print() const	//prints the system (coefficients only)
{
    for (int i = 0; i < height; i++) {
        cout << "| ";
        for (int j = 0; j < width - 1; j++) {
            cout << fixed << setw(10) << coefficients[i][j] << " ";
        }
        cout << "| " << fixed << setw(10) << coefficients[i][width - 1] << " |" << endl;
    }
    cout << endl;
}

void Matrix::prn_solutions(bool sci)	 //prints the solutions only
{
    int err = 0;
    if (solutions == nullptr)
        err = solve();
    if (err == 0) {
        for (int i = 0; i < height; i++) {
            cout << (sci ? scientific : fixed) << "solutions" << i + 1 << " = " << solutions[i] << ";" << endl;
        }
    }
    else {
        error(err);
    }
}

void Matrix::error(int err) const	// prints an error based on coefficients code
{
    switch (err) {
        case 1:
            cout << "The row diverges!" << endl;
            break;
        case 2:
            cout << "The system cannot be solved because of zeros!" << endl;
            break;
        default:
            cout << "Unknown error" << endl;
            break;
    }
}

int main()
{
    char file[] = "system.txt";
    Matrix* m = new Matrix(file);
    m->print();
    m->prn_solutions(false);
}
