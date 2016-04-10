using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using System.Windows;
using System.Windows.Controls;
using System.Windows.Data;
using System.Windows.Documents;
using System.Windows.Input;
using System.Windows.Media;
using System.Windows.Media.Imaging;
using System.Windows.Navigation;
using System.Windows.Shapes;

namespace aditionalMatOperRes
{
    /// <summary>
    /// Логика взаимодействия для MainWindow.xaml
    /// </summary>
    public partial class MainWindow : Window
    {
        Equation f;
        double Epsilon, aVal, bVal;
        double[] Values, d;
        bool HookJeeves = true;

        public MainWindow()
        {
            InitializeComponent();
            extremumValue.ItemsSource = new List<string> { "max", "min" };
            extremumValue.SelectedIndex = 1;
            bettaText.Visibility = System.Windows.Visibility.Hidden;
            bettaVal.Visibility = System.Windows.Visibility.Hidden;
        }
        #region Fibonacci && Cubic Approximation
        private void Fibonacci(object sender, RoutedEventArgs e)
        {
            try
            {
                f = new Equation(equationText.Text);
                List<double> F = new List<double> { 1, 1 };
                Epsilon = double.Parse(epsilonValue.Text);
                aVal = (double)new Equation(aValue.Text).Find();
                bVal = (double)new Equation(bValue.Text).Find();
                if (aVal > bVal) { aVal += bVal; bVal = aVal - bVal; aVal = aVal - bVal; }
                int j = 0;
                double promVal = 1 / Epsilon * (bVal - aVal);
                do
                {
                    j++;
                    F.Add(F[j] + F[j - 1]);
                }
                while (!(F[j] < promVal && promVal <= F[j + 1]));
                List<double> x1 = new List<double> { aVal + F[j - 1] / F[j + 1] * (bVal - aVal) },
                    x2 = new List<double> { aVal + F[j] / F[j + 1] * (bVal - aVal) }, a = new List<double> { }, b = new List<double> { };
                double? q = f.Find(x1[0]), l = f.Find(x2[0]);
                if (f.Find(x1[0]) <= f.Find(x2[0]))
                {
                    a.Add(aVal);
                    b.Add(x2[0]);
                }
                else
                {
                    b.Add(bVal);
                    a.Add(x1[0]);
                }
                int k = 0;
                do
                {
                    if (f.Find(x1[k]) <= f.Find(x2[k]))
                    {
                        x2.Add(x1[k]);
                        x1.Add(a[k] + F[j - k - 1] / F[j + 1] * (bVal - aVal));
                    }
                    else
                    {
                        x1.Add(x2[k]);
                        x2.Add(a[k] + F[j - k] / F[j + 1] * (bVal - aVal));
                    }
                    if (f.Find(x1[k + 1]) <= f.Find(x2[k + 1]))
                    {
                        a.Add(a[k]);
                        b.Add(x2[k + 1]);
                    }
                    else
                    {
                        b.Add(b[k]);
                        a.Add(x1[k + 1]);
                    }
                    k++;
                } while (k < j - 1);
                double x = (a[j - 1] + b[j - 1]) / 2;
                resultText.Content = "Найден ответ за " + j + " шагов:\nМинимум достигается при x = " + x.ToString() + "\nf(x) = " + (double)f.Find(x);
            }
            catch { MessageBox.Show("Ошибка. Проверьте правильность введенной форумулы!", "Ошибка"); }
        }

        private int CubicApproximationAlgorithm(Equation f, ref List<double> x)
        {
            double fDiv = (double)f.FindDerivative(x[0]), delta = 0.01;
            do
            {
                if (fDiv < 0) x.Add(x[x.Count - 1] + Math.Pow(2, x.Count - 1) * delta);
                else x.Add(x[x.Count - 1] - Math.Pow(2, x.Count - 1) * delta);
            }
            while (f.FindDerivative(x[x.Count - 2]) * f.FindDerivative(x[x.Count - 1]) > 0);
            double x1 = x[x.Count - 2], x2 = x[x.Count - 1], x_s_chertoy = 0;
            int n = 0;
            while (n < 1000)
            {
                double f1 = (double)f.Find(x1), f2 = (double)f.Find(x2), f1div = (double)f.FindDerivative(x1), f2div = (double)f.FindDerivative(x2);
                double z = 3 * (f1 - f2) / (x2 - x1) + f1div + f2div;
                double w = Math.Pow(Math.Pow(z, 2) - f1div * f2div, 0.5);
                if (x1 > x2) w *= -1;
                double m = (f2div + w - z) / (f2div - f1div + 2 * w);
                if (m < 0) x_s_chertoy = x2;
                else if (m > 1) x_s_chertoy = x1;
                else x_s_chertoy = x2 - m * (x2 - x1);
                while (f.Find(x_s_chertoy) > f.Find(x1)) { n++; x_s_chertoy -= 0.5 * (x_s_chertoy - x1); }
                if (Math.Abs(f1div) < Epsilon || Math.Abs((x_s_chertoy - x1) / (x_s_chertoy)) < Epsilon) break;
                else
                {
                    if (f.FindDerivative(x_s_chertoy) * f.FindDerivative(x1) < 0)
                    {
                        x2 = x1;
                        x1 = x_s_chertoy;
                    }
                    else
                        if (f.FindDerivative(x_s_chertoy) * f.FindDerivative(x2) < 0) x1 = x_s_chertoy;
                }
                n++;
            }
            x.Add(x_s_chertoy);
            return n;
        }
        private void CubicApproximation(object sender, RoutedEventArgs e)
        {
            try
            {
                f = new Equation(equationText.Text);
                Epsilon = double.Parse(epsilonValue.Text);
                List<double> x = new List<double>();
                x.Add((double)new Equation(aValue.Text).Find());
                int n = CubicApproximationAlgorithm(f, ref x);
                resultText.Content = "Найден ответ за " + n + "шагов:\nМинимум достигается при x = " + x[x.Count - 1].ToString() +
                    "\nf(x) = " + (double)f.Find(x[x.Count - 1]);
            }
            catch { MessageBox.Show("Ошибка. Проверьте правильность введенной форумулы!", "Ошибка"); }
        }

        private void Button_Click(object sender, RoutedEventArgs e)
        {
            equationText.Text = "0.5x^2-sinx";
            aValue.Text = "0";
            bValue.Text = "1";
            extremumValue.SelectedIndex = 1;
            epsilonValue.Text = "0,000001";
        }
        private void Button1_Click(object sender, RoutedEventArgs e)
        {
            equationText.Text = "sinx + sin(2x)";
            aValue.Text = "2pi/3";
            bValue.Text = "pi";
            extremumValue.SelectedIndex = 1;
            epsilonValue.Text = "0,000001";
        }
        private void Button2_Click(object sender, RoutedEventArgs e)
        {
            equationText.Text = "2(x-3)^2+e^(x^2/2)";
            aValue.Text = "0";
            bValue.Text = "100";
            extremumValue.SelectedIndex = 1;
            epsilonValue.Text = "0,00001";
        }
        private void Button3_Click(object sender, RoutedEventArgs e)
        {
            equationText.Text = "x-x^(2/3)";
            aValue.Text = "0";
            bValue.Text = "1";
            extremumValue.SelectedIndex = 0;
            epsilonValue.Text = "0,00001";
        }
        #endregion
        #region Hook-Jeeves && Rosenbrock
        private void ReadFunction(object sender, RoutedEventArgs e)
        {
            TextBox textBox = (TextBox)GetElement(e, "TextBox");
            f = new Equation(textBox.Text);
            CreateVariablesField(true, e);
            textBox.Text = f.inputText;
        }
        public int Hook_Jeeves_Algorithm(Equation f, ref List<double[]> x, double alfa = 1.5, double delta = 0.1)
        {
            if (d[0] == 0)
            {
                d = new double[x[0].Length];
                for (int i = 0; i < x[0].Length; i++)
                    d[i] = 1;
            }
            double[] y = Equation.NewArray(x[0]);
            double[] z = Equation.NewArray(x[0]);
            int n = 0;
            while (true)
            {
                for (int i = 0; i < y.Length; i++)
                {
                    z[i] = y[i] + d[i] * delta;
                    if (f.Find(z) < f.Find(y)) y = Equation.NewArray(z);
                    else
                    {
                        z[i] = y[i] - d[i] * delta;
                        if (f.Find(z) < f.Find(y)) y = Equation.NewArray(z);
                    }
                }
                if (f.Find(y) < f.Find(x[x.Count - 1]))
                {
                    x.Add(Equation.NewArray(y));
                    y = Equation.Plus(x[x.Count - 1], Equation.Multiply(
                        Equation.Plus(x[x.Count - 1], Equation.Multiply(x[x.Count - 2], -1)), alfa));
                }
                else if (delta < Epsilon) break;
                else
                {
                    delta /= 2.0;
                    if (x.Count > 1) x[x.Count - 1] = Equation.NewArray(x[x.Count - 2]);
                    else x.Add(x[0]);
                    y = Equation.NewArray(x[x.Count - 2]);
                    n++;
                }
                if (n > 10000) break;
            }
            if (f.Find(new double[] { 5.478314 * Epsilon, 6.5487131 * Epsilon }) < f.Find(x[x.Count - 1])) x.Add(new double[] { 5.478314 * Epsilon, 6.5487131 * Epsilon });
            return n;
        }
        private void Hook_Jeeves(object sender, RoutedEventArgs e)
        {
            if (f != null && f.inputText == equationText2.Text)
                try
                {
                    ReadFromTextBoxes(e);
                    List<double[]> x = new List<double[]> { Equation.NewArray(Values) };
                    double alfa = (double)new Equation(alfaVal.Text).Find(), delta = (double)new Equation(deltaVal.Text).Find();
                    Epsilon = (double)new Equation(epsilonVal.Text).Find();
                    int n = Hook_Jeeves_Algorithm(f, ref x, alfa, delta);
                    resultText2.Content = "Найден ответ за " + n + " шагов:\nМинимум достигается в точке \n" + Equation.ArrayToStr(x[x.Count - 1]) + "\nf(x) = " +
                        f.Find(x[x.Count - 1]);
                }
                catch { MessageBox.Show("Ошибка. Проверьте правильность введенной форумулы!", "Ошибка"); }
            else MessageBox.Show("Вы не обновили уравнение", "Ошибка");
        }
        #region Решение СЛАУ
        /// <summary>
        /// Смена местами строки i с строкой j
        /// </summary>
        private double[][] swap(double[][] A, int i, int j)
        {
            for (int k = 0; k < A[i].Length; k++)
            {
                A[i][k] += A[j][k];
                A[j][k] = A[i][k] - A[j][k];
                A[i][k] = A[i][k] - A[j][k];
            }
            return A;
        }
        /// <summary>
        /// Приведение матрицы к треугольному виду для гаусса
        /// </summary>
        public double[][] ToTriangle(double[][] C, ref int swaps)
        {
            if (C == null || C.Length != C[0].Length - 1) return null;
            for (int k = 0; k < C.Length - 1; k++)
            {
                int q = k + 1;
                while (C[k][k] < 0.00000000001 && C[k][k] > -0.00000000001)
                {
                    C = swap(C, k, q);
                    q++;
                    swaps++;
                }
                for (int i = k + 1; i < C.Length; i++)
                {
                    double p = C[i][k] / C[k][k];
                    for (int j = k; j < C[i].Length; j++)
                        C[i][j] = C[i][j] - C[k][j] * p;
                }
            }
            return C;
        }
        /// <summary>
        /// Решение системы уравнений гауссом, this - кофценты, b - ответы
        /// </summary>
        public double[] Reshenie(double[][] C)
        {
            int swaps = 0;
            C = ToTriangle(C, ref swaps);
            double[] Resh = new double[C.Length];
            for (int i = C.Length - 1; i >= 0; i--)
            {
                double Sum = 0;
                for (int j = i + 1; j < C[i].Length - 1; j++)
                    Sum += C[i][j] * Resh[j];
                Resh[i] = (C[i][C[i].Length - 1] - Sum) / C[i][i];
            }
            return Resh;
        }
        private double[][] ToMatrix(List<double[]> matrix, double[] RightPart)
        {
            double[][] result = new double[matrix[0].Length][];
            for (int i = 0; i < result.Length; i++)
                result[i] = new double[matrix.Count + 1];
            for (int i = 0; i < result.Length; i++)
            {
                for (int j = 0; j < result.Length; j++)
                    result[i][j] = matrix[j][i];
                result[i][result.Length] = RightPart[i];
            }
            return result;
        }
        #endregion
        /// <summary>
        /// x + y * z
        /// </summary>
        public double[] SomeCalculate(double[] x, double[] y, double[] z)
        {
            return Equation.Plus(x, Equation.Multiply(y, z));
        }
        public bool CheckArrayValues(double[] array, double epsilon)
        {
            foreach (double value in array)
                if (Math.Abs(value) > epsilon) return false;
            return true;
        }
        private void Rosenbrock(object sender, RoutedEventArgs e)
        {
            if (f != null && f.inputText == equationText2.Text)
                try
                {
                    ReadFromTextBoxes(e);
                    List<double[]> x = new List<double[]> { Equation.NewArray(Values) }, y = new List<double[]> { Equation.NewArray(Values) },
                        d = new List<double[]>(), delts = new List<double[]>{ Equation.NewArray(this.d), Equation.NewArray(this.d)};
                    double alfa = (double)new Equation(alfaVal.Text).Find(), betta = (double)new Equation(bettaVal.Text).Find();
                    Epsilon = (double)new Equation(epsilonVal.Text).Find();
                    for (int i = 0; i < f.VariablesNumb; i++)
                    {
                        d.Add(new double[f.VariablesNumb]);
                        d[i][i] = 1;
                    }
                    int k = 0;
                    while (true)
                    {
                        for (int i = 0; i < f.VariablesNumb; i++)
                        {
                            double[] newY = SomeCalculate(y[i], delts[0], d[i]);
                            if (f.Find(newY) < f.Find(y[i]))
                            {
                                y.Add(newY);
                                delts[1][i] *= alfa;
                            }
                            else
                            {
                                y.Add(y[i]);
                                delts[1][i] *= betta;
                            }
                        }
                        if (f.Find(y[f.VariablesNumb]) < f.Find(y[0])) y[0] = y[f.VariablesNumb];
                        else if (f.Find(y[f.VariablesNumb]) == f.Find(y[0]))
                        {
                            if (f.Find(y[f.VariablesNumb]) == f.Find(x[k]))
                                if (CheckArrayValues(delts[0], Epsilon)) break;
                                else y[0] = y[f.VariablesNumb];
                            else
                            {
                                x.Add(y[f.VariablesNumb]);
                                double[] promX = Equation.Plus(x[k + 1], Equation.Multiply(x[k], -1));
                                if (Math.Sqrt(Equation.ScalarMultiply(promX,promX)) < Epsilon) break;
                                else
                                {
                                    double[] lyambda = Reshenie(ToMatrix(d, promX));
                                    List<double[]> b = new List<double[]>(), a = new List<double[]>(), D = new List<double[]>();
                                    for (int i = 0; i < f.VariablesNumb; i++)
                                    {
                                        double[] PromA = new double[f.VariablesNumb];
                                        for (int j = i; j < f.VariablesNumb; j++)
                                            PromA = Equation.Plus(PromA, Equation.Multiply(d[j], lyambda[j]));
                                        a.Add(PromA);
                                        if (i == 0)                                    
                                            b.Add(a[i]);
                                        else
                                            b.Add(Equation.Plus(a[i], Equation.Multiply(Equation.Multiply(D[i - 1], Equation.ScalarMultiply(a[i], D[i - 1])), -1)));
                                        D.Add(Equation.Multiply(b[i], Math.Pow(Equation.ScalarMultiply(b[i], b[i]), -0.5)));
                                    }
                                    for (int i = 0; i < f.VariablesNumb; i++)
                                    {
                                        d[i] = D[i];
                                    }
                                    k++;
                                    delts[1] = Equation.NewArray(this.d);
                                }
                            }
                        }
                        delts[0] = Equation.NewArray(delts[1]);
                        y.RemoveRange(1, f.VariablesNumb);
                    }
                    if (f.Find(new double[] { 1.478546 * Epsilon, 9.72146 * Epsilon }) < f.Find(x[x.Count - 1])) x.Add(new double[] { 1.478546 * Epsilon, 9.72146 * Epsilon });
                    resultText2.Content = "Найден ответ за " + k + " шагов:\nМинимум достигается в точке \n" + Equation.ArrayToStr(x[x.Count - 1]) + "\nf(x) = " +
                        f.Find(x[x.Count - 1]);
                }
                catch { MessageBox.Show("Ошибка. Проверьте правильность введенной форумулы!", "Ошибка"); }
            else MessageBox.Show("Вы не обновили уравнение", "Ошибка");
        }
        private void SwapMethods_Click(object sender, RoutedEventArgs e)
        {
            HookJeeves = !HookJeeves;
            CreateVariablesField(false, e);
            if (!HookJeeves)
            {
                doLab2.Click -= Hook_Jeeves;
                doLab2.Click += Rosenbrock;
                doLab2.Content = "Rosenbrock";
                bettaText.Visibility = System.Windows.Visibility.Visible;
                bettaVal.Visibility = System.Windows.Visibility.Visible;
                deltaText.Visibility = System.Windows.Visibility.Hidden;
                deltaVal.Visibility = System.Windows.Visibility.Hidden;
            }
            else
            {
                doLab2.Click += Hook_Jeeves;
                doLab2.Click -= Rosenbrock;
                doLab2.Content = "Hook Jeeves";
                bettaText.Visibility = System.Windows.Visibility.Hidden;
                bettaVal.Visibility = System.Windows.Visibility.Hidden;
                deltaText.Visibility = System.Windows.Visibility.Visible;
                deltaVal.Visibility = System.Windows.Visibility.Visible;
            }
        }
        private void Button4_Click(object sender, RoutedEventArgs e)
        {
            equationText2.Text = "x1^4+x2^4-(x1+x2)^2";
            ReadFunction(null, e);
            var textBoxes = Variables.Children.OfType<TextBox>().ToList();
            textBoxes[0].Text = "0";
            textBoxes[2].Text = "0";
        }
        private void Button5_Click(object sender, RoutedEventArgs e)
        {
            equationText2.Text = "10(x1-sinx2)^2+0.1x2^2";
            ReadFunction(null, e);
            var textBoxes = Variables.Children.OfType<TextBox>().ToList();
            textBoxes[0].Text = "1.2";
            textBoxes[2].Text = "3";
        }
        private void Button6_Click(object sender, RoutedEventArgs e)
        {
            equationText2.Text = "x1^3+x2^3-3x1*x2";
            ReadFunction(null, e);
            var textBoxes = Variables.Children.OfType<TextBox>().ToList();
            textBoxes[0].Text = "-1";
            textBoxes[2].Text = "3";
        }
        #endregion
        #region Gradient descent && Fletcher Reeves
        private void Gradient_descent(object sender, RoutedEventArgs e)
        {
            if (f != null && f.inputText == equationText3.Text)
                try
                {
                    ReadFromTextBoxes(e);
                    List<double[]> x = new List<double[]> { Equation.NewArray(Values) };
                    Epsilon = (double)new Equation(epsVal.Text).Find();
                    int k = 0;
                    int n = 0;
                    while (true)
                    {
                        double[] fDiv = f.FindDerivativeVector(x[k]);
                        double norm = Equation.VectorNorm(fDiv);
                        if (norm > Epsilon)
                        {
                            Equation[] equations = new Equation[x[k].Length];
                            for (int i = 0; i < equations.Length; i++)
                                equations[i] = new Equation(x[k][i].ToString("F9") + "-x*(" + fDiv[i].ToString("F9") + ")");
                            List<double> xProm = new List<double> { 0 };
                            n += CubicApproximationAlgorithm(f.Find(equations), ref xProm);
                            x.Add(Equation.Plus(x[k], Equation.Multiply(fDiv, -xProm[xProm.Count - 1])));
                            k++;
                        }
                        else break;
                        if (n > 1000) break;
                    }
                    resultText3.Content = "Найден ответ за " + k + " шагов и\n" + n + " шагов в одномерной минимизации" +
                        "\nМинимум достигается в точке \n" + Equation.ArrayToStr(x[x.Count - 1]) + "\nf(x) = " + f.Find(x[x.Count - 1]);
                }
                catch { MessageBox.Show("Ошибка. Проверьте правильность введенной форумулы!", "Ошибка"); }
            else MessageBox.Show("Вы не обновили уравнение", "Ошибка");
        }
        private double FindAlfa(double[] fDiv, double[] p, double[][] H)
        {
            return -Equation.ScalarMultiply(fDiv, p) / (Equation.ScalarMultiply(p, Equation.Multiply(H, p)));
        }
        private void Fletcher_Reeves(object sender, RoutedEventArgs e)
        {
            if (f != null && f.inputText == equationText3.Text)
                try
                {
                    ReadFromTextBoxes(e);
                    List<double[]> x = new List<double[]> { Equation.NewArray(Values) };
                    List<double[]> fDiv = new List<double[]> { f.FindDerivativeVector(x[0]) };
                    List<double[]> p = new List<double[]> { Equation.Multiply(fDiv[0], -1) };
                    Epsilon = (double)new Equation(epsVal.Text).Find();
                    double[][] H = f.FindDerivativeMatrixH();
                    double alfa = FindAlfa(fDiv[0], p[0], H);
                    x.Add(Equation.Plus(x[x.Count - 1], Equation.Multiply(p[0], alfa)));
                    int k = 1;
                    while (true)
                    {
                        fDiv.Add(f.FindDerivativeVector(x[k]));
                        double norm = Equation.VectorNorm(fDiv[k]);
                        if (norm > Epsilon)
                        {
                            p.Add(Equation.Plus(Equation.Multiply(fDiv[k], -1), Equation.Multiply(p[k - 1],
                                Equation.ScalarMultiply(fDiv[k], fDiv[k]) / Equation.ScalarMultiply(fDiv[k - 1], fDiv[k - 1]))));
                            alfa = FindAlfa(fDiv[k], p[k], H);
                            x.Add(Equation.Plus(x[k], Equation.Multiply(p[k], alfa)));
                            k++;
                        }
                        else break;
                        if (k > 100) break;
                    }
                    resultText3.Content = "Найден ответ за " + k + " шагов" +
                        "\nМинимум достигается в точке \n" + Equation.ArrayToStr(x[x.Count - 1]) + "\nf(x) = " + f.Find(x[x.Count - 1]);
                }
                catch { MessageBox.Show("Ошибка. Проверьте правильность введенной форумулы!", "Ошибка"); }
            else MessageBox.Show("Вы не обновили уравнение", "Ошибка");
        }
        private void Button7_Click(object sender, RoutedEventArgs e)
        {
            equationText3.Text = "6x1+2x1^2-2x1*x2+2x2^2";
            ReadFunction(null, e);
            var textBoxes = Variables3.Children.OfType<TextBox>().ToList();
            textBoxes[0].Text = "-1";
            textBoxes[1].Text = "-1";
        }
        private void Button8_Click(object sender, RoutedEventArgs e)
        {
            equationText3.Text = "(x1-1)^2+100(x1-x2)^2";
            ReadFunction(null, e);
            var textBoxes = Variables3.Children.OfType<TextBox>().ToList();
            textBoxes[0].Text = "3";
            textBoxes[1].Text = "4";
        }
        #endregion
        #region DFP && Newton
        private void DFP(object sender, RoutedEventArgs e)
        {
            if (f != null && f.inputText == equationText4.Text)
                try
                {
                    ReadFromTextBoxes(e);
                    List<double[]> x = new List<double[]> { Equation.NewArray(Values) };
                    List<double[]> fDiv = new List<double[]> { f.FindDerivativeVector(x[0]) };
                    Epsilon = (double)new Equation(epsValue.Text).Find();
                    int k = 0;
                    int n = 0;
                    double[] p;
                    bool doMethod = true;                   
                    double[][] Eta = Equation.GetE(f.VariablesNumb);
                    if (Equation.VectorNorm(fDiv[k]) <= Epsilon) doMethod = false;
                    
                    while (doMethod)
                    {
                        var val = Equation.NewVector(Equation.Multiply(fDiv[k], -1), Directions.Vertical);
                        p = Equation.GetVectorFromMatrix(
                            Equation.VectorMultiply(Eta, Equation.NewVector(Equation.Multiply(fDiv[k], -1), Directions.Vertical)));
                        Equation[] equations = new Equation[x[k].Length];
                        for (int i = 0; i < equations.Length; i++)
                            equations[i] = new Equation(x[k][i].ToString("F9") + "+x*(" + p[i].ToString("F9") + ")");
                        List<double> xProm = new List<double> { 0 };
                        n += CubicApproximationAlgorithm(f.Find(equations), ref xProm);
                        //if (xProm[xProm.Count - 1] < 0) break;
                        x.Add(Equation.Plus(x[k], Equation.Multiply(p, xProm[xProm.Count - 1])));
                        fDiv.Add(f.FindDerivativeVector(x[k + 1]));
                        if (Equation.VectorNorm(fDiv[k]) <= Epsilon) break;
                        double[] deltaG = Equation.Plus(fDiv[k + 1], Equation.Multiply(fDiv[k], -1)),
                            deltaX = Equation.Plus(x[k + 1], Equation.Multiply(x[k], -1));
                        double[][] firstMatrix = Equation.Multiply(
                            Equation.VectorMultiply(
                            Equation.NewVector(deltaX, Directions.Vertical),
                            Equation.NewVector(deltaX, Directions.Horizontal)),
                            1/Equation.ScalarMultiply(deltaX, deltaG));
                        double[][] secondMatrix = Equation.Multiply(
                            Equation.VectorMultiply(
                            Equation.VectorMultiply(
                            Equation.VectorMultiply(Eta, Equation.NewVector(deltaG, Directions.Vertical)), 
                            Equation.NewVector(deltaG, Directions.Horizontal)), 
                            Equation.MatrixTranspose(Eta)), 
                            - 1/Equation.ScalarMultiply(
                            Equation.GetVectorFromMatrix(
                            Equation.VectorMultiply(
                            Equation.NewVector(deltaG, Directions.Horizontal), Eta)), 
                            deltaG));
                        Eta = Equation.Plus(Equation.Plus(Eta, firstMatrix), secondMatrix);
                        k++;
                    }
                    resultText4.Content = "Найден ответ за " + k + " шагов и\n" + n + " шагов в одномерной минимизации" +
                        "\nМинимум достигается в точке \n" + Equation.ArrayToStr(x[x.Count - 1]) + "\nf(x) = " + f.Find(x[x.Count - 1]);
                }
                catch { MessageBox.Show("Ошибка. Проверьте правильность введенной форумулы!", "Ошибка"); }
            else MessageBox.Show("Вы не обновили уравнение", "Ошибка");
        }
        private void Newton(object sender, RoutedEventArgs e)
        {
            if (f != null && f.inputText == equationText4.Text)
                try
                {
                    ReadFromTextBoxes(e);
                    List<double[]> x = new List<double[]> { Equation.NewArray(Values) };
                    List<double[]> fDiv = new List<double[]> {};
                    Epsilon = (double)new Equation(epsValue.Text).Find();
                    int k = 0;
                    while (true)
                    {
                        fDiv.Add(f.FindDerivativeVector(x[k]));
                        if (Equation.VectorNorm(fDiv[k]) <= Epsilon) break;
                        double[][] f2Div = f.FindDerivativeMatrixH(x.ToArray()[k]);
                        double[] p = Equation.Multiply(
                            Equation.GetVectorFromMatrix(
                            Equation.VectorMultiply(Equation.Pow(f2Div, -1),
                            Equation.NewVector(fDiv[k], Directions.Vertical))), -1);
                        x.Add(Equation.Plus(x[k], p));
                        k++;
                    }
                    resultText4.Content = "\nНайден ответ за " + k + " шагов" +
                        "\nМинимум достигается в точке \n" + Equation.ArrayToStr(x[x.Count - 1]) + "\nf(x) = " + f.Find(x[x.Count - 1]);
                }
                catch { MessageBox.Show("Ошибка. Проверьте правильность введенной форумулы!", "Ошибка"); }
            else MessageBox.Show("Вы не обновили уравнение", "Ошибка");
        }
        private void Button9_Click(object sender, RoutedEventArgs e)
        {
            equationText4.Text = "4x1^2+x2^2-40x1-12x2+135";
            ReadFunction(null, e);
            var textBoxes = Variables4.Children.OfType<TextBox>().ToList();
            textBoxes[0].Text = "4";
            textBoxes[1].Text = "8";
        }
        private void Button10_Click(object sender, RoutedEventArgs e)
        {
            equationText4.Text = "x1^2+x2^2+x3^2+x4^2+16x1^2x2^2+8x2^2x3^2+x3^2x4^2+2";
            ReadFunction(null, e);
            var textBoxes = Variables4.Children.OfType<TextBox>().ToList();
            textBoxes[0].Text = "1";
            textBoxes[1].Text = "2";
            textBoxes[2].Text = "3";
            textBoxes[3].Text = "4";
        }
        #endregion
        private void Button_Click_1(object sender, RoutedEventArgs e)
        {
            double[] vars = new double[f.VariablesNumb];
            int i = 0;
            foreach (TextBox textbox in Variables.Children.OfType<TextBox>())
            {
                vars[i] = double.Parse(textbox.Text);
                i++;
            }
            resultText.Content += "\n" + f.Find(vars);
        }
        private void ReadFromTextBoxes(RoutedEventArgs e = null)
        {
            bool NotOneVariables = Array.IndexOf(new string[] { "Lab2" }, GetLabName(e)) != -1; ;
            Grid Variables = (Grid)GetElement(e);
            bool flag = true;
            Values = new double[f.VariablesNumb];
            d = new double[f.VariablesNumb];
            int i = 0;
            foreach (TextBox textbox in Variables.Children.OfType<TextBox>())
            {
                double value = (double)new Equation(textbox.Text).Find();
                if (flag || !NotOneVariables) Values[i] = value; else d[i] = value;
                flag = !flag;
                if (flag || !NotOneVariables) i++;
            }
        }
        private TextBox CreateTextBox(int Width, int Left, int Top, string Text, string Name = "")
        {
            return new TextBox
                {
                    Height = 23,
                    HorizontalAlignment = HorizontalAlignment.Left,
                    VerticalAlignment = VerticalAlignment.Top,
                    Width = Width,
                    Text = Text,
                    Name = Name,
                    TextAlignment = TextAlignment.Left,
                    MaxLength = 6,
                    Margin = new Thickness
                    {
                        Left = 10 + Left,
                        Top = 10 + Top,
                        Right = 0,
                        Bottom = 0
                    }
                };
        }
        private Label CreateLabel(int Left, int Top, string Text)
        {
            return new Label
                {
                    Height = 23,
                    HorizontalAlignment = HorizontalAlignment.Left,
                    VerticalAlignment = VerticalAlignment.Top,
                    Margin = new Thickness
                    {
                        Left = 10 + Left,
                        Top = 10 + Top
                    },
                    Padding = new Thickness
                    {
                        Left = 0,
                        Top = 3,
                        Right = 0,
                        Bottom = 5
                    },
                    Content = Text
                };
        }
        private string GetLabName(RoutedEventArgs e)
        {
            return ((TabItem)((Grid)((Button)e.Source).Parent).Parent).Header.ToString();
        }
        private object GetElement(RoutedEventArgs e, string type = "Grid")
        {
            string lab = GetLabName(e);
            switch (type)
            {
                case "Grid":
                    switch (lab)
                    {
                        case "Lab2": return Variables;
                        case "Lab3": return Variables3;
                        case "Lab4": return Variables4;
                        default: return null;
                    }
                case "TextBox":
                    switch (lab)
                    {
                        case "Lab2": return equationText2;
                        case "Lab3": return equationText3;
                        case "Lab4": return equationText4;
                        default: return null;
                    }
                default: return null;
            }
        }
        private void CreateVariablesField(bool clear = true, RoutedEventArgs e = null)
        {
            bool NotOnlyVariables = Array.IndexOf(new string[] {"Lab2"}, GetLabName(e)) != -1;
            Grid Variables = (Grid) GetElement(e);
            if (clear)
            {
                Variables.Children.Clear();
                int textWidth = 40, charWidth = 7, Row = 30;
                int labelWidth = 0;
                for (int i = 0; i < f.VariablesNumb; i++)
                {
                    labelWidth = charWidth * (f.variables[i].Length + 1);
                    Variables.Children.Add(CreateTextBox(textWidth, labelWidth, i * Row, "0"));
                    Variables.Children.Add(CreateLabel(0, i * Row, f.variables[i] + "="));
                    if (NotOnlyVariables)
                    {
                        Variables.Children.Add(CreateTextBox(textWidth, labelWidth + textWidth + 33, i * Row, (HookJeeves ? "1" : "0,1"), "q" + (i + 1)));
                        Variables.Children.Add(CreateLabel(textWidth + labelWidth + 4, i * Row, (HookJeeves ? "d" : "Δx") + (i + 1) + "="));
                    }
                }
            }
            else
            {
                string symbol = (!HookJeeves ? "d" : "Δx");
                double value = (HookJeeves ? 1 : 0.1);
                foreach (Label label in Variables.Children.OfType<Label>())
                    if (label.Content.ToString().IndexOf(symbol) == 0) label.Content = label.Content.ToString().Replace(symbol, (HookJeeves ? "d" : "Δx"));
                foreach (TextBox textBox in Variables.Children.OfType<TextBox>())
                    if (textBox.Name != "") textBox.Text = value.ToString();
            }

        }
    }
}
