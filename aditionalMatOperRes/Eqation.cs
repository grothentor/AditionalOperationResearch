using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace aditionalMatOperRes
{
    public enum Directions { Horizontal, Vertical}
    public class Equation
    {
        public string inputText { set; get; }
        public int VariablesNumb { set; get; }
        public List<string> variables { set; get; }
        private List<object> RPN_Text { set; get; }

        public static Dictionary<char, int> Signs = new Dictionary<char, int> 
        { { '+', 2 }, { '-', 2 }, { '*', 4 }, { '/', 4 }, { '^', 6 }, { '(', 0 }, { ')', 0 } };
        public static int FunctionsPriority = 3;
        public static List<char> Numbers = new List<char> { '1', '2', '3', '4', '5', '6', '7', '8', '9', '0', ',' };
        public static List<string> Functions = new List<string> { "sin", "cos", "tg", "ctg" };
        double dx = 0.00001, dx2 = 0.000005;
        public Equation()
        {
            inputText = "";
            RPN_Text = new List<object>();
            VariablesNumb = 0;
            variables = new List<string>();
        }
        public Equation(string Equation, int variant = 1)
        {
            Parse(Equation);
            FindVariables();
        }
        public void SetRPN(List<object> newRPN)
        {
            this.RPN_Text = new List<object>(newRPN.ToArray());
        }
        public static int Check_Priority(char a, object b)
        {
            if (b is string) return Signs[a] - FunctionsPriority;
            if (Signs.Keys.Contains(a) && Signs.Keys.Contains((char)b))
                return Signs[a] - Signs[(char)b];
            else return 2;
        }
        public static string Check_Function(string text)
        {
            foreach (string function in Functions)
                if (text.Contains(function)) return function;
            return null;
        }
        public void RebuildInputText()
        {
            string result = null;
            foreach (char val in inputText)
                if (val != ' ' && (int)val != 32) if (result != null) result += val; else result = val.ToString();
            result = result.Replace("(-", "(0-");
            result = result.Replace(".", ",");
            if (result[0] == '-') result = "0-" + result.Remove(0, 1);
            inputText = result;
        }
        public int CheckFunction(string text)
        {
            for (int i = text.Length - 1; i >= 0; i--)
                if (Functions.Contains(text.Substring(i))) return i;
            return -1;
        }
        public bool Parse(string Equation)
        {
            inputText = Equation;
            RebuildInputText();
            List<object> RPN_Text = new List<object>();
            List<object> Stack = new List<object>();
            string currentNumber = "";
            string currentVariable = "";
            foreach (char Symbol in inputText)
            {
                int value = CheckFunction(currentVariable);
                if (value >= 0)
                {
                    Stack.Add(currentVariable.Remove(0, value));
                    currentVariable = currentVariable.Remove(value, currentVariable.Length - value);
                    if (currentVariable != "")
                    {
                        RPN_Text.Add(currentVariable);
                        currentVariable = "";
                        while (Stack.Count != 0 && Check_Priority('*', Stack[Stack.Count - 1]) <= 0)
                        {
                            RPN_Text.Add(Stack[Stack.Count - 1]);
                            Stack.RemoveAt(Stack.Count - 1);
                        }
                        Stack.Add('*');
                    }                    
                }
                if (Numbers.Contains(Symbol) && currentVariable == "") currentNumber += Symbol;
                else
                {
                    if (Signs.Keys.Contains(Symbol))
                    {                        
                        if (currentNumber != "")
                        {
                            RPN_Text.Add(double.Parse(currentNumber));
                            currentNumber = "";
                            if (Symbol == '(')
                            {
                                while (Stack.Count != 0 && Check_Priority('*', Stack[Stack.Count - 1]) <= 0)
                                {
                                    RPN_Text.Add(Stack[Stack.Count - 1]);
                                    Stack.RemoveAt(Stack.Count - 1);
                                }
                                Stack.Add('*');
                            }
                        }
                        if (currentVariable != "")
                        {
                            RPN_Text.Add(currentVariable);
                            currentVariable = "";
                            if (Symbol == '(')
                            {
                                while (Stack.Count != 0 && Check_Priority('*', Stack[Stack.Count - 1]) <= 0)
                                {
                                    RPN_Text.Add(Stack[Stack.Count - 1]);
                                    Stack.RemoveAt(Stack.Count - 1);
                                }
                                Stack.Add('*');
                            }
                        }
                        if (Symbol == ')')
                        {
                            while (!(Stack[Stack.Count - 1] is char) || (char)Stack[Stack.Count - 1] != '(')
                            {
                                RPN_Text.Add(Stack[Stack.Count - 1]);
                                Stack.RemoveAt(Stack.Count - 1);
                            }
                            Stack.RemoveAt(Stack.Count - 1);
                        }
                        else
                        {                            
                            while (Symbol != '(' && Stack.Count != 0 && Check_Priority(Symbol, Stack[Stack.Count - 1]) <= 0)
                            {
                                RPN_Text.Add(Stack[Stack.Count - 1]);
                                Stack.RemoveAt(Stack.Count - 1);
                            }
                            Stack.Add(Symbol);
                        }
                    }
                    else
                    {
                        if (currentNumber != "")
                        {
                            RPN_Text.Add(double.Parse(currentNumber));
                            currentNumber = "";
                            while (Stack.Count != 0 && Check_Priority('*', Stack[Stack.Count - 1]) <= 0)
                            {
                                RPN_Text.Add(Stack[Stack.Count - 1]);
                                Stack.RemoveAt(Stack.Count - 1);

                            }
                            Stack.Add('*');
                        }                        
                        currentVariable += Symbol;
                    }
                }
            }
            if (currentNumber != "")
            {
                RPN_Text.Add(double.Parse(currentNumber));
                currentNumber = "";
            }
            if (currentVariable != "")
            {
                RPN_Text.Add(currentVariable);
                currentVariable = "";
            }
            while (Stack.Count != 0)
            {
                RPN_Text.Add(Stack[Stack.Count - 1]);
                Stack.RemoveAt(Stack.Count - 1);
            }
            this.RPN_Text = RPN_Text;
            return true;
        }
        private void FindVariables()
        {
            variables = new List<string>();
            for (int i = 0; i < RPN_Text.Count; i++)
                if (RPN_Text[i] is string && !Functions.Contains(RPN_Text[i]) && !variables.Contains(RPN_Text[i]))
                    if (RPN_Text[i].ToString() == "e") RPN_Text[i] = Math.E;
                    else if (RPN_Text[i].ToString() == "pi") RPN_Text[i] = Math.PI;
                    else variables.Add((string)RPN_Text[i]);
            VariablesNumb = variables.Count;
            variables.Sort();
        }
        public override string ToString()
        {
            string result = "";
            foreach (object val in RPN_Text)
                result += val.ToString();
            return result;
        }
        #region ArrayOperations
        public static double[] Plus( double[] array1, double[] array2)
        {
            double[] newArray = new double[array1.Length];
            if (array1.Length != array2.Length) return null;
            for (int i = 0; i < array1.Length; i++)
                newArray[i] += array1[i] + array2[i];
            return newArray;
        }
        public static double[][] Plus(double[][] array1, double[][] array2)
        {
            double[][] newArray = new double[array1.Length][];
            if (array1.Length != array2.Length) return null;
            for (int i = 0; i < array1.Length; i++)
                newArray[i] = Plus(array1[i], array2[i]);
            return newArray;
        }
        public static double[] Plus(double[] array1, double delta)
        {
            double[] newArray = new double[array1.Length];
            for (int i = 0; i < array1.Length; i++)
                 newArray[i] = array1[i] + delta;
            return newArray;
        }
        public static double[] Multiply(double[] array1, double delta)
        {
            double[] newArray = new double[array1.Length];
            for (int i = 0; i < array1.Length; i++)
                newArray[i] = array1[i] * delta;
            return newArray;
        }
        public static double[][] Multiply(double[][] array1, double delta)
        {
            double[][] newArray = new double[array1.Length][];
            for (int i = 0; i < array1.Length; i++)
                newArray[i] = Multiply(array1[i], delta);
            return newArray;
        }
        public static double[] Multiply(double[][] array1, double[] array2)
        {
            double[] result = new double[array1.Length];
            for (int i = 0; i < array1.Length; i++)
                result[i] = ScalarMultiply(array1[i], array2);
            return result;
        }
        public static double ScalarMultiply(double[] array1, double[] array2)
        {
            double Summ = 0;
            for (int i = 0; i < array1.Length; i++)
                Summ += array1[i] * array2[i];
            return Summ;
        }
        public static double VectorNorm(double[] array)
        {
            double Summ = 0;
            for (int i = 0; i < array.Length; i++)
                Summ += Math.Pow(array[i], 2);
            return Math.Sqrt(Summ);
        }
        public static double[] Pow(double[] array1, double power)
        {
            double[] newArray = new double[array1.Length];
            for (int i = 0; i < array1.Length; i++)
                newArray[i] = Math.Pow(array1[i], power);
            return newArray;
        }
        public static double[] Multiply(double[] array1, double[] array2)
        {
            double[] newArray = new double[array1.Length];
            if (array1.Length != array2.Length) return null;
            for (int i = 0; i < array1.Length; i++)
                newArray[i]= array1[i] * array2[i];
            return newArray;
        }
        public static double[][] VectorMultiply(double[][] array1, double[][] array2)
        {
            double[][] newArray = new double[array1.Length][];
            for (int i = 0; i < array1.Length; i++)
                newArray[i] = new double[array2[0].Length];
            if (array1[0].Length != array2.Length) return null;
            for (int i = 0; i < array1.Length; i++)
                for (int j = 0; j < array1[i].Length; j++)
                    for (int k = 0; k < array2[j].Length; k++)
                        newArray[i][k] += array1[i][j]*array2[j][k];
            return newArray;
        }
        public static double[] NewArray(double[] array)
        {
            double[] newArray = new double[array.Length];
            for (int i = 0; i < array.Length; i++) newArray[i] = array[i];
            return newArray;
        }
        public static double[][] NewMatrix(double[][] array)
        {
            double[][] newArray = new double[array.Length][];
            for (int i = 0; i < array.Length; i++) newArray[i] = NewArray(array[i]);
            return newArray;
        }
        public static double[][] NewVector(double[] array, Directions direction)
        {
            double[][] newVector = new double[][] { };
            switch (direction)
            {
                case Directions.Horizontal:
                    newVector = new double[1][];
                    newVector[0] = NewArray(array);
                    break;
                case Directions.Vertical:
                    newVector = new double[array.Length][];
                    for (int i = 0; i < array.Length; i++)
                        newVector[i] = new double[] { array[i] };
                    break;
            }
            return newVector;
        }
        public static double[] GetVectorFromMatrix(double[][] matrix)
        {
            double[] array = new double[]{};
            if (matrix.Length == 1)
                array = NewArray(matrix[0]);
            else
                if (matrix.Length < 1 || matrix[0].Length != 1) return null;
                else
                {
                    array = new double[matrix.Length];
                    for (int i = 0; i < matrix.Length; i++)
                        array[i] = matrix[i][0];
                }
            return array;
        }
        public static double[][] MatrixTranspose(double[][] matrix)
        {
            if (matrix.Length == 0 || matrix.Length != matrix[0].Length) return null;
            double[][] newMatrix = new double[matrix.Length][];
            for (int i = 0; i < matrix.Length; i++)
            {
                newMatrix[i] = new double[matrix.Length];
                for (int j = 0; j < matrix.Length; j++)
                    newMatrix[i][j] = matrix[j][i];
            }
            return newMatrix;
        }
        public static string ArrayToStr(double[] array)
        {
            string result = "(";
            for (int i = 0; i < array.Length; i++)
                result += array[i].ToString("F6") + ";";
            result = result.Remove(result.Length - 1);
            result += ")";
            return result;
        }
        public static string MatrixToStr(double[][] array)
        {
            string result = "";
            for (int i = 0; i < array.Length; i++)
                result += ArrayToStr(array[i]) + "\n";
            return result;
        }
        public static double[][] GetE(int n)
        {
            double[][] matrix = new double[n][];
            for (int i = 0; i < n; i++)
            {
                matrix[i] = new double[n];
                matrix[i][i] = 1;
            }
            return matrix;
        }
        public static double[][] Pow(double[][] matrix, int n)
        {
            if (n < -1) return null;
            if (n == -1) return Multiply(MatrixTranspose(MatrixAddition(matrix)), (double)(1 / MatrixDeterminant(matrix)));
            if (n == 0) return GetE(matrix.Length);
            double[][] result = VectorMultiply(matrix, matrix);
            for (int i = 2; i < n; i++)
            {
                result = VectorMultiply(result, matrix);
            }
            return result;
        }
        public static double MatrixDeterminant(double[][] matrix){
            double Sum = 1;
            int swaps = 0;
            double[][] C = MatrixToTriangle(ref swaps, matrix);
            for (int i = 0; i < matrix.Length; i++)
                Sum *= C[i][i];
            return swaps % 2 == 0 ? Sum : -Sum;
        }
        public static double[][] MatrixToTriangle(ref int swaps, double[][] matrix)
        {
            double[][] C = NewMatrix(matrix);
            for (int k = 0; k < C.Length - 1; k++)
            {
                int q = k + 1;
                while (C[k][k] < 0.00000000001 && C[k][k] > -0.00000000001)
                {
                    swap(C, k, q);
                    q++;
                    swaps++;
                }
                for (int i = k + 1; i < C.Length; i++)
                {
                    double p = C[i][k] / C[k][k];
                    for (int j = k; j < C.Length; j++)
                        C[i][j] = C[i][j] - C[k][j] * p;
                }
            }
            return C;
        }
        private static void swap(double[][] matrix, int i, int j)
        {
            for (int k = 0; k < matrix.Length; k++)
            {
                matrix[i][k] += matrix[j][k];
                matrix[j][k] = matrix[i][k] - matrix[j][k];
                matrix[i][k] = matrix[i][k] - matrix[j][k];
            }
        }
        public static double[][] MatrixAddition(double[][] matrix)
        {
            double[][] C = NewMatrix(matrix);
            for (int i = 0; i < matrix.Length; i++)
                for (int j = 0; j < matrix[i].Length; j++)
                {
                    double[][] a = MatrixDelete(matrix, i, j);
                    double det = (double)MatrixDeterminant(a);
                    C[i][j] = Math.Pow(-1, i + j) * det;
                }
            return C;
        }
        public static double[][] MatrixDelete(double[][] matrix, int Row, int Col)
        {
            double[][] C = new double[matrix.Length-1][];
            int k = 0, n = 0;
            for (int i = 0; i < C.Length; i++)
            {
                C[i] = new double[matrix[i].Length - 1];
                if (i < Row) k = i; else k = i + 1;
                for (int j = 0; j < C[i].Length; j++)
                {
                    if (j < Col) n = j; else n = j + 1;
                    C[i][j] = matrix[k][n];
                }
            }
            return C;
        }
        #endregion
        public double[] FindDerivativeVector(params double[] variables)
        {
            decimal dx = (decimal)this.dx;
            double[] Deveretives = new double[variables.Length];
            for (int i = 0; i < variables.Length; i++)
            {
                double[] variables2 = (double[])variables.Clone();
                variables2[i] = variables[i] + this.dx;
                Deveretives[i] = (double)(((decimal)Find(variables2) - (decimal)Find(variables)) / dx);
            }
            return Deveretives;
        }
        public double FindSecondDerivative(double[] delta1, double[] delta2, double[] x = null)
        {
            decimal dx = (decimal)this.dx;
            decimal dx2 = (decimal)this.dx2;
            if ( x == null) x = new double[variables.Count];
            decimal f1 = (decimal)Find(Plus(x, Plus(delta1, delta2))), 
                f2 = (decimal)Find(Plus(x, Plus(delta2, Multiply(delta1, -1)))),
                f3 = (decimal)Find(Plus(x, Plus(delta1, Multiply(delta2, -1)))),
                f4 = (decimal)Find(Plus(x, Multiply(Plus(delta2, delta1), -1)));
            decimal f1Div = (f1 - f2) / (2 * dx),
                f2Div = (f3 - f4) / (2 * dx),
                result = (f1Div - f2Div) / (2 * dx2);
            return (double)result;
        }
        public double[][] FindDerivativeMatrixH(double[] x = null)
        {
            double[][] Deveretives = new double[variables.Count][];
            double[] delta = new double[variables.Count];
            for (int i = 0; i < variables.Count; i++)
            {
                Deveretives[i] = new double[variables.Count];
                for (int j = 0; j < variables.Count; j++)
                {
                    double[] delta1 = (double[])delta.Clone();
                    double[] delta2 = (double[])delta.Clone();
                    delta1[i] += dx;
                    delta2[j] += dx2;
                    Deveretives[i][j] = FindSecondDerivative(delta1, delta2, x);
                }
            }
            return Deveretives;
        }
        public double? FindDerivative(params double[] variables)
        {
            double[] variables2 = new double[variables.Length];
            for (int i = 0; i < variables.Length; i++)
                variables2[i] = variables[i] + dx;
            return (Find(variables2) - Find(variables)) / dx;
        }
        public double? FindMultiDerivative(int order = 1, params double[] variables)
        {
            if (order < 0) return null;
            if (order == 0) return Find(variables);
            if (order == 1) return FindDerivative(variables);
            double[] variables2 = new double[variables.Length];
            for (int i = 0; i < variables.Length; i++)
                variables2[i] = variables[i] + dx;
            return (FindMultiDerivative(order - 1, variables2) - FindMultiDerivative(order - 1, variables)) / dx;
        }
        public double? Find(params double[] variables)
        {
            if (variables.Length != VariablesNumb) return null;
            List<object> RPN_Text = new List<object>(this.RPN_Text.ToArray());
            for (int i = 0; i < RPN_Text.Count; i++)
                if (RPN_Text[i] is string && !Functions.Contains(RPN_Text[i])) 
                    RPN_Text[i] = variables[this.variables.IndexOf((string)RPN_Text[i])];
            while (RPN_Text.Count != 1)
            {
                for (int i = 0; i < RPN_Text.Count; i++)
                    if (RPN_Text[i] is char)
                    {
                        DoOperation(ref RPN_Text, i);
                        break;
                    } 
                    else if (RPN_Text[i] is string)
                    {
                        DoFunction(ref RPN_Text, i);
                        break;
                    }
            }
            return (double)RPN_Text[0];
        }
        public Equation Find(Equation[] equations)
        {
            if (equations.Length != VariablesNumb) return null;
            List<object> RPN_Text = new List<object>(this.RPN_Text.ToArray());
            for (int i = 0; i < RPN_Text.Count; i++)
                if (RPN_Text[i] is string && !Functions.Contains(RPN_Text[i]))
                {
                    List<object> newRPN = equations[this.variables.IndexOf((string)RPN_Text[i])].RPN_Text;
                    RPN_Text.RemoveAt(i);
                    RPN_Text.InsertRange(i, newRPN);
                    i += newRPN.Count - 1;
                }
            Equation equation = new Equation();
            equation.SetRPN(RPN_Text);
            equation.FindVariables();
            return equation;
        }
        private void DoFunction(ref List<object> RPN_Text, int n)
        {
            switch ((string)RPN_Text[n])
            {
                case "sin": RPN_Text[n - 1] = Math.Sin((double)RPN_Text[n - 1]); break;
                case "cos": RPN_Text[n - 1] = Math.Cos((double)RPN_Text[n - 1]); break;
                case "tg": RPN_Text[n - 1] = Math.Tan((double)RPN_Text[n - 1]); break;
                case "ctg": RPN_Text[n - 1] = 1/Math.Tan((double)RPN_Text[n - 1]); break;
            }
            RPN_Text.RemoveAt(n);
        }
        private void DoOperation(ref List<object> RPN_Text, int n)
        {
            switch ((char)RPN_Text[n])
            {
                case '+': RPN_Text[n - 2] = (double)RPN_Text[n - 2] + (double)RPN_Text[n - 1]; break;
                case '-': RPN_Text[n - 2] = (double)RPN_Text[n - 2] - (double)RPN_Text[n - 1]; break;
                case '*': RPN_Text[n - 2] = (double)RPN_Text[n - 2] * (double)RPN_Text[n - 1]; break;
                case '/': RPN_Text[n - 2] = (double)RPN_Text[n - 2] / (double)RPN_Text[n - 1]; break;
                case '^': RPN_Text[n - 2] = Math.Pow((double)RPN_Text[n - 2], (double)RPN_Text[n - 1]); break;
            }
            RPN_Text.RemoveRange(n - 1, 2);
        }
    }
}
