namespace Lab6
{
    public class Blue
    {
        public delegate int Finder(int[,] matrix, out int row, out int col);
        public delegate void SortRowsStyle(int[,] matrix, int row);
        public delegate void ReplaceMaxElements(int[,] matrix, int row, int maxValue);
        public delegate int[] GetTriangle(int[,] matrix);
        public void Task1(ref int[,] matrix)
        {

            // code here
            if (matrix.GetLength(0) == matrix.GetLength(1))
            {
                int index = FindDiagonalMaxIndex(matrix);
                RemoveRow(ref matrix, index);
            }
            // end

        }
        public int Task2(int[,] A, int[,] B, int[,] C)
        {
            int answer = 0; // 1 - increasing   0 - no sequence   -1 - decreasing

            // code here
            double[] matrix = new double[3];
            matrix[0] = GetAverageExceptEdges(A);
            matrix[1] = GetAverageExceptEdges(B);
            matrix[2] = GetAverageExceptEdges(C);
            if (matrix[0] > matrix[1] && matrix[1] > matrix[2]) answer = -1;
            if (matrix[0]<matrix[1] && matrix[1]<matrix[2]) answer = 1;
            // end

            return answer;
        }
        public void Task3(ref int[,] matrix, Func<int[,], int> method)
        {

            // code here
            if (matrix.GetLength(0) == matrix.GetLength(1))
            {
                int colindex = method(matrix);
                RemoveColumn(ref matrix, colindex);
            }
            // end

        }
        public void Task4(ref int[,] matrix)
        {

            // code here
            for (int j = matrix.GetLength(1) - 1; j >= 0; j--)
            {
                if (CheckZerosInColumn(matrix, j) == false) RemoveColumn(ref matrix, j);
            }
            // end

        }
        public void Task5(ref int[,] matrix, Finder find)
        {

            // code here
            int row, col;
            int chislo = find(matrix, out  row, out  col);
            while (find(matrix, out  row, out  col) == chislo) RemoveRow(ref matrix, row);
            // end

        }
        public void Task6(int[,] matrix, SortRowsStyle sort)
        {

            // code here
            for (int i = 0; i < matrix.GetLength(0); i++)
            {
                if (i%3==0)
                {
                    sort(matrix, i);
                }
            }
            // end

        }
        public void Task7(int[,] matrix, ReplaceMaxElements transform)
        {

            // code here
            for (int i = 0; i < matrix.GetLength(0); i++)
            {
                int max = FindMaxInRow(matrix, i);
                transform(matrix, i, max);
            }
            // end

        }
        public double[,] Task8(double a, double b, double h, Func<double, double> getSum, Func<double, double> getY)
        {
            double[,] answer = null;

            // code here
            answer = GetSumAndY(a, b, h, getSum, getY);
            // end

            return answer;
        }
        public int Task9(int[,] matrix, GetTriangle triangle)
        {
            int answer = 0;

            // code here
            if (matrix.GetLength(0) == matrix.GetLength(1))
            {
                int[] array = triangle(matrix);
                answer = Sum(array);
            }
            // end

            return answer;
        }
        public bool Task10(int[][] array, Predicate<int[][]> func)
        {
            bool res = false;

            // code here
            res = func(array);  
            // end
            return res;
        }
        
        public int FindDiagonalMaxIndex(int[,] matrix)
        {
            int max = matrix[0, 0];
            int index = 0;
            for (int i = 0; i < matrix.GetLength(0); i++)
            {
                if (matrix[i, i] > max)
                {
                    max = matrix[i, i];
                    index = i;
                }
            }
            return index;
        }

        public void RemoveRow(ref int[,] matrix, int rowIndex)
        {
            int[,] temp = new int[matrix.GetLength(0)-1, matrix.GetLength(1)];
            for (int i = 0; i < rowIndex; i++)
            {
                for (int j = 0; j < matrix.GetLength(1); j++)
                {
                    temp[i, j] = matrix[i, j];
                }
            }
            for (int i = rowIndex+1; i < matrix.GetLength(0); i++)
            {
                for (int j = 0; j < matrix.GetLength(1); j++)
                {
                    temp[i-1, j] = matrix[i, j];
                }
            }
            matrix = new int[temp.GetLength(0), temp.GetLength(1)];
            matrix = temp;
        }

        public double GetAverageExceptEdges(int[,] matrix)
        {
            double avg = 0;
            int max, row, col;
            int min, mrow, mcol;
            max = FindMax(matrix, out row, out col);
            min = FindMin(matrix, out mrow, out mcol);
            for (int i=0;i<matrix.GetLength(0);i++)
            {
                for (int j = 0; j < matrix.GetLength(1); j++)
                {
                    avg += matrix[i, j];
                }
            }

            avg -= (max + min);
            avg /= matrix.Length-2;
            return avg;
        }

        public int FindUpperColIndex(int[,] matrix)
        {
            int max = matrix[0, 1];
            int index = 0;
            for (int i=0;i<matrix.GetLength(0);i++)
            {
                for (int j = i+1; j < matrix.GetLength(1); j++)
                {
                    if (matrix[i, j] > max)
                    {
                        max = matrix[i, j];
                        index = j;
                    }
                }
            }
            return index;
        }

        public int FindLowerColIndex(int[,] matrix)
        {
            int max = matrix[0, 0];
            int index = 0;
            for (int j=0;j<matrix.GetLength(1);j++)
            {
                for (int i = j; i < matrix.GetLength(0); i++)
                {
                    if (matrix[i, j] > max)
                    {
                        max = matrix[i, j];
                        index = j;
                    }
                }
            }
            return index;
        }

        public void RemoveColumn(ref int[,] matrix, int col)
        {
            int[,] temp = new int[matrix.GetLength(0), matrix.GetLength(1)-1];
            for (int j = 0; j < col; j++)
            {
                for (int i = 0; i < matrix.GetLength(0); i++)
                {
                    temp[i, j] = matrix[i, j];
                }
            }
            for (int j = col+1; j < matrix.GetLength(1); j++)
            {
                for (int i = 0; i < matrix.GetLength(0); i++)
                {
                    temp[i, j-1] = matrix[i, j];
                }
            }
            matrix = new int[temp.GetLength(0), temp.GetLength(1)];
            matrix = temp;
        }

        public bool CheckZerosInColumn(int[,] matrix, int col)
        {
            bool res = false;
            for (int i = 0; i < matrix.GetLength(0); i++)
            {
                if (matrix[i, col] == 0) res = true;
            }
            return res;
        }
        
        public int FindMax(int[,] matrix, out int row, out int col)
        {
            int max = int.MinValue;
            row = 0;
            col = 0;
            for (int i=0;i<matrix.GetLength(0);i++)
            {
                for (int j = 0; j < matrix.GetLength(1); j++)
                {
                    if (matrix[i, j] > max)
                    {
                        max = matrix[i, j];
                        row = i;
                        col = j;
                    }
                }
            }
            return max;
        }

        public int FindMin(int[,] matrix, out int row, out int col)
        {
            int min = int.MaxValue;
            row = 0;
            col = 0;
            for (int i=0;i<matrix.GetLength(0);i++)
            {
                for (int j = 0; j < matrix.GetLength(1); j++)
                {
                    if (matrix[i, j] < min)
                    {
                        min = matrix[i, j];
                        row = i;
                        col = j;
                    }
                }
            }
            return min;
        }

        public void SortRowAscending(int[,] matrix, int row)
        {
            for (int i = 0; i < matrix.GetLength(1); i++)
            {
                for (int j = 0; j < matrix.GetLength(1)-1-i; j++)
                {
                    if (matrix[row, j] > matrix[row, j + 1])
                    {
                        (matrix[row, j], matrix[row, j + 1]) = (matrix[row, j+1], matrix[row, j]);
                    }
                }
            }
        }

        public void SortRowDescending(int[,] matrix, int row)
        {
            for (int i = 0; i < matrix.GetLength(1); i++)
            {
                for (int j = 0; j < matrix.GetLength(1)-1-i; j++)
                {
                    if (matrix[row, j] < matrix[row, j + 1])
                    {
                        (matrix[row, j], matrix[row, j + 1]) = (matrix[row, j+1], matrix[row, j]);
                    }
                }
            }
        }

        public int FindMaxInRow(int[,] matrix, int row)
        {
            int max = int.MinValue;
            for (int j = 0; j < matrix.GetLength(1); j++)
            {
                if (matrix[row, j] > max) max = matrix[row, j];
            }
            return max;
        }

        public void ReplaceByZero(int[,] matrix, int row, int maxValue)
        {
            for (int j = 0; j < matrix.GetLength(1); j++)
            {
                if (matrix[row, j] == maxValue) matrix[row, j] = 0;
            }
        }

        public void MultiplyByColumn(int[,] matrix, int row, int maxValue)
        {
            for (int j = 0; j < matrix.GetLength(1); j++)
            {
                if (matrix[row, j] == maxValue) matrix[row, j] *= j+1;
            }
        }

        public double[,] GetSumAndY(double a, double b, double h, Func<double, double> sum, Func<double, double> y)
        {
            double[,] matrix = new double[(int)((b - a)/h)+1, 2];
            for (int i = 0; i < matrix.GetLength(0); i++)
            {
                double x = a + i * h;
                matrix[i, 0] = sum(x);
                matrix[i, 1] = y(x);
            }
            return matrix;
        }

        public double SumA(double x)
        {
            double sum = 1.0;        
            double prevSum = 0.0;
            double factorial = 1.0;
    
            int i = 1;
            do
            {
                prevSum = sum;
                factorial *= i;
                sum += Math.Cos(i * x) / factorial;
                i++;
            }
            while (Math.Abs(sum - prevSum) >= 0.0001);  
    
            return sum;
        }

        public double YA(double x)
        {
            double res = 0;
            res = Math.Exp(Math.Cos(x))*Math.Cos(Math.Sin(x));
            return res;
        }

        public double SumB(double x)
        {
            double s = -2.0 * Math.PI * Math.PI / 3.0;

            for (double i = 1; ; i += 1.0)
            {
                s += Math.Pow(-1, i) * Math.Cos(i * x) / (i * i);

                if (Math.Abs(Math.Pow(-1, i) * Math.Cos(i * x) / (i * i)) < 0.000001) break;
            }
            return s;
        }

        public double YB(double x)
        {
            return (x * x) / 4.0 - 3.0 * (Math.PI * Math.PI) / 4.0;
        }

        public int Sum(int[] array)
        {
            int sum = 0;
            for (int i = 0; i < array.Length; i++)
            {
                sum += array[i]*array[i];
            }
            return sum;
        }

        public int[] GetSum(GetTriangle transformer, int[,] matrix)
        {
            int[] res = transformer(matrix);
            return res;
        }

        public int[] GetUpperTriangle(int[,] matrix)
        {
            int n = matrix.GetLength(0);
            int[] array = new int[n*(n+1)/2];
            int k = 0;
            for (int i = 0; i < matrix.GetLength(0); i++)
            {
                for (int j = i; j < matrix.GetLength(1); j++)
                {
                    array[k++] = matrix[i, j];
                }
            }
            return array;;
        }
        public int[] GetLowerTriangle(int[,] matrix)
        {
            int n = matrix.GetLength(0);
            int[] array = new int[n*(n+1)/2];
            int k = 0;
            for (int j = 0; j < matrix.GetLength(1); j++)
            {
                for (int i = j; i < matrix.GetLength(0); i++)
                {
                    array[k++] = matrix[i, j];
                }
            }
            return array;
        }

        public bool CheckTransformAbility(int[][] array)
        {
            bool res = false;
            int count = 0;
            for (int i = 0; i < array.Length; i++)
            {
                count += array[i].Length;
            }

            if (count % array.Length == 0) res = true;
            return res;
        }

        public bool CheckSumOrder(int[][] array)
        {
            if (array.Length < 2) return true;
            bool res = false;
            int[] temp = new int[array.Length];
            int k = 0;
            for (int i = 0; i < array.Length; i++)
            {
                for (int j = 0; j < array[i].Length; j++)
                {
                    temp[k] += array[i][j];
                }
                k++;
            }

            for (int i = 0; i < temp.Length-1; i++)
            {
                if (temp[i] > temp[i + 1])
                {
                    res = true;
                }
                else
                {
                    res = false;
                    break;
                }
                
            }
            if (res) return res;
            for (int i = 0; i < temp.Length - 1; i++)
            {
                if (temp[i] < temp[i + 1])
                {
                    res = true;
                }
                else
                {
                    res = false;
                    break;
                }
            }
            return res;
        }

        public bool CheckArraysOrder(int[][] array)
        {
            
            bool res = false;
            for (int i = 0; i < array.Length; i++)
            {
                 if (array[i].Length < 2) return true;
                for (int j = 0; j < array[i].Length-1; j++)
                {
                    if (array[i][j] >= array[i][j+1])
                    {
                        res = true;
                    }
                    else
                    {
                        res = false;
                        break;
                    }
                    
                }
                if (res) return true;
                for (int j = 0; j < array[i].Length-1; j++)
                {
                    if (array[i][j] <= array[i][j+1])
                    {
                        res = true;
                    }
                    else
                    {
                        res = false;
                        break;
                    }
                    
                }
                if (res) return true;
            }
            return res;
        }
       
    }
}
