#if 0
typedef struct Complex
{
	double R;
	double C;
}Complex;

Complex* add(Complex* S1, const Complex* S2) {
	S1->R = S1->R + S2->R;
	S1->C = S1->C + S2->R;
	return S1;
}
Complex* subtraction(Complex* S1, const Complex* S2) {
	S1->R = S1->R - S2->R;
	S1->C = S1->C - S2->C;
	return S1;
}
Complex* multiplicaion(Complex* S1, const Complex* S2) {
	S1->R = S1->R * S2->R - S1->C * S2->C;
	S1->C = S1->C * S2->R + S2->C * S1->R;
	return S1;
}
Complex* division(Complex* S1, Complex* S2) {
	double ZN = (S1->R * S1->R) + (S1->C * S2->C);

	S2->C = (S2->C) * (-1);
	multiplicaion(S1, S2);
	(S1->R) / ZN; (S1->C) / ZN;
	return S1;
}
int printC(const Complex* S1)
{
	printf("Real part: %lf\n", S1->R);
	printf("Imaginary part: %lf\n", S1->C);
	return 0;
}

int readC(const Complex* num) {
	printf("Enter real part(dot and numbers from 0 to 9): ");
	scanf("%lf", &num->R);
	printf("Enter imaginary part(dot and numbers from 0 to 9): ");
	scanf("%lf", &num->C);
	return 0;
}

void complex_calc() {

	Complex S1, S2;
	S1.R = 0; S1.C = 0; S2.R = 0; S2.C = 0;

	static char L = 'E';
	while (1) {

		printf("%s", "Enter Action('-','+','/','*','E') <");

		if (0 == scanf(" %c", &L))abort();
		getchar();
		if (L == 'E')
		{
			return 0;
		}
		readC(&S1);
		readC(&S2);
		switch (L)
		{
		case '+': printC(add(&S1, &S2)); break;
		case '-': printC(subtraction(&S1, &S2)); break;
		case '*': printC(multiplicaion(&S1, &S2)); break;
		case '/': printC(division(&S1, &S2)); break;
		default:abort();
		}

		
	}
	printf("%s", "ïîñëåäíåå ÷èñëî");
	printC(&S1);







	system("pause");

}
#endif