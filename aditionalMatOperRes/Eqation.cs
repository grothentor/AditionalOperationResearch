using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace aditionalMatOperRes
{     
    public class Equation
    {
        public string inputText { set; get; }
        public int VariablesNumb { set; get; }
        public List<string> variables { set; get; }
        public List<object> RPN_Text { set; get; }
        public static Dictionary<char, int> Signs = new Dictionary<char, int> 
        { { '+', 2 }, { '-', 2 }, { '*', 4 }, { '/', 4 }, { '^', 6 }, { '(', 0 }, { ')', 0 } };
        public static int FunctionsPriority = 3;
        public static List<char> Numbers = new List<char> { '1', '2', '3', '4', '5', '6', '7', '8', '9', '0', ',' };
        public static List<string> Functions = new List<string> { "sin", "cos", "tg", "ctg" };
        double dx = 0.000005;
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
        public static double ScalarMultiply(double[] array1, double[] array2)
        {
            double Summ = 0;
            for (int i = 0; i < array1.Length; i++)
                Summ += array1[i] * array2[i];
            return Summ;
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
        public static double[] VectorMultiply(double[] array1, double[] array2)
        {
            double[] newArray = new double[array1.Length];
            if (array1.Length != array2.Length) return null;
            for (int i = 0; i < array1.Length; i++)
                for (int j = 0; j < array1.Length; j++)
                newArray[i] += array1[i] * array2[i];
            return newArray;
        }
        public static double[] NewArray(double[] array)
        {
            double[] newArray = new double[array.Length];
            for (int i = 0; i < array.Length; i++) newArray[i] = array[i];
            return newArray;
        }
        public static string ArrayToStr(double[] array)
        {
            string result = "(";
            for (int i = 0; i < array.Length; i++)
                result += array[i] + ";";
            result = result.Remove(result.Length - 1);
            result += ")";
            return result;
        }
        #endregion
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
