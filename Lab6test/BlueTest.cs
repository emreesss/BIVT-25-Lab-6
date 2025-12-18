using System.Transactions;

namespace Lab6test
{
    [TestClass]
    public sealed class BlueTest
    {
        Lab6.Blue _main = new Lab6.Blue();
        const double E = 0.0001;
        Data _data = new Data();

        [TestMethod]
        public void Test01()
        {
            // Arrange
            var input = _data.GetMatrixes();
            var answer = new int[][,] {
                new int[,] {
                    {1, 2, 3, 4, 5, 6, 7},
                    {5, 6, 7, 8, 9, 10, 11},
                    {9, 10, 11, 12, 13, 14, 15},
                    {13, 14, 15, 16, 17, 18, 19},
                    {0, 1, 2, 3, 4, 5, 6},
                },
                new int[,] {
                    {1},
                    {5},
                    {9},
                    {13},
                },
                new int[,] {
                    {1, 2, 3, 4, 5, 6},
                    {5, 6, 7, 8, 9, 11},
                    {0, 2, 3, 4, 5, 6},
                },
                new int[,] {
                    {1, 2, 4, 6},
                    {5, -6, 7, 11},
                    {-1, 4, -5, 6},
                },
                new int[,] {
                    {1, 2, 3, 4, 5, 6, 7},
                    {5, 6, 7, 8, -9, 10, 11},
                    {9, 10, -11, -12, -13, -14, -15},
                    {-13, -14, 15, 16, 17, 18, -19},
                },
                new int[,] {
                    {1, 2, 3},
                    {0, -2, -3},
                },
                new int[,] {
                    {-9, -10, -11, -14, -15, 6},
                },
                new int[,] {
                    {1, 2, 3, 4, -5, -6, -7},
                    {5, 11, -17, 11, -10, 6, 5},
                    {-9, -10, -11, -14, -15, -16, 1},
                    {-9, -10, -11, -14, -15, -6, -2},
                    {-9, -10, -11, -14, -15, 6, 4},
                },
                new int[,] {
                    {1, 2, 3, 4, -5, -6, -7},
                    {5, 11, -17, 11, -10, 6, 5},
                    {-9, -10, -11, -14, -15, -16, 1},
                    {-9, -10, -11, -14, 15, -6, -2},
                    {-9, -10, -11, -14, -15, 6, 4},
                    {5, 11, -17, 11, -10, 6, -5},
                    {1, 1, -2, 3, -4, 0, 0},
                    {0, -2, -3, -4, -5, 0, 5},
                },
                new int[,] {
                    {1, 2, 3, 4, -5},
                    {-9, -10, -11, -14, -15},
                    {-9, -10, -11, -14, -6},
                    {0, -2, -3, -4, -5},
                }
            };
            // Act
            for (int i = 0; i < answer.Length; i++)
            {
                _main.Task1(ref input[i]);
            }
            // Assert
            for (int i = 0; i < answer.Length; i++)
            {
                Assert.AreEqual(answer[i].GetLength(0), input[i].GetLength(0));
                for (int j = 0; j < answer[i].GetLength(0); j++)
                {
                    Assert.AreEqual(answer[i].GetLength(1), input[i].GetLength(1));
                    for (int k = 0; k < answer[i].GetLength(1); k++)
                    {
                        Assert.AreEqual(answer[i][j, k], input[i][j, k]);
                    }
                }
            }
        }
        [TestMethod]
        public void Test02()
        {
            // Arrange
            var inputA = _data.GetMatrixes();
            var inputB = _data.GetMatrixes();
            var inputC = _data.GetMatrixes();
            Array.Copy(inputA, 0, inputB, 1, 9);
            Array.Reverse(inputC);
            Array.Copy(inputC, 0, inputB, 1, 3);
            var answer = new int[] { 0, 0, -1, -1, 0, 0, 1, 0, 0, 1 };
            var test = new int[answer.Length];
            // Act
            for (int i = 0; i < answer.Length; i++)
            {
                test[i] = _main.Task2(inputA[i], inputB[i], inputC[i]);
            }
            // Assert
            Assert.AreEqual(answer.Length, test.Length);
            for (int i = 0; i < answer.Length; i++)
            {
                Assert.AreEqual(answer[i], test[i]);
            }
        }
        [TestMethod]
        public void Test03a()
        {
            // Arrange
            var input = _data.GetMatrixes();
            var answer = new int[][,] {
                new int[,] {
                    {1, 2, 3, 4, 5, 6, 7},
                    {5, 6, 7, 8, 9, 10, 11},
                    {9, 10, 11, 12, 13, 14, 15},
                    {13, 14, 15, 16, 17, 18, 19},
                    {0, 1, 2, 3, 4, 5, 6},
                },
                new int[,] {
                    {1},
                    {5},
                    {9},
                    {13},
                },
                new int[,] {
                    {1, 2, 3, 4, 5, 6},
                    {5, 6, 7, 8, 9, 11},
                    {0, 2, 3, 4, 5, 6},
                },
                new int[,] {
                    {1, 2, 4},
                    {5, -6, 7},
                    {-1, 4, -5},
                    {1, 4, 5},
                },
                new int[,] {
                    {1, 2, 3, 4, 5, 6, 7},
                    {5, 6, 7, 8, -9, 10, 11},
                    {9, 10, -11, -12, -13, -14, -15},
                    {-13, -14, 15, 16, 17, 18, -19},
                },
                new int[,] {
                    {1, 2},
                    {5, 11},
                    {0, -2},
                },
                new int[,] {
                    {-9, -10, -11, -14, -15, 6},
                },
                new int[,] {
                    {1, 2, 3, 4, -5, -6, -7},
                    {5, 11, -17, 11, -10, 6, 5},
                    {-9, -10, -11, -14, -15, -16, 1},
                    {-9, -10, -11, -14, -15, -6, -2},
                    {-9, -10, -11, -14, -15, 6, 4},
                },
                new int[,] {
                    {1, 2, 3, 4, -5, -6, -7},
                    {5, 11, -17, 11, -10, 6, 5},
                    {-9, -10, -11, -14, -15, -16, 1},
                    {-9, -10, -11, -14, 15, -6, -2},
                    {-9, -10, -11, -14, -15, 6, 4},
                    {5, 11, -17, 11, -10, 6, -5},
                    {1, 1, -2, 3, -4, 0, 0},
                    {0, -2, -3, -4, -5, 0, 5},
                },
                new int[,] {
                    {1, 2, 3, -5},
                    {5, 11, -17, 7},
                    {-9, -10, -11, -15},
                    {-9, -10, -11, -6},
                    {0, -2, -3, -5},
                }
            };
            // Act
            for (int i = 0; i < answer.Length; i++)
            {
                _main.Task3(ref input[i], _main.FindUpperColIndex);
            }
            // Assert
            for (int i = 0; i < answer.Length; i++)
            {
                Assert.AreEqual(answer[i].GetLength(0), input[i].GetLength(0));
                for (int j = 0; j < answer[i].GetLength(0); j++)
                {
                    Assert.AreEqual(answer[i].GetLength(1), input[i].GetLength(1));
                    for (int k = 0; k < answer[i].GetLength(1); k++)
                    {
                        Assert.AreEqual(answer[i][j, k], input[i][j, k]);
                    }
                }
            }
        }
        [TestMethod]
        public void Test03b()
        {
            // Arrange
            var input = _data.GetMatrixes();
            var answer = new int[][,] {
                new int[,] {
                    {1, 2, 3, 4, 5, 6, 7},
                    {5, 6, 7, 8, 9, 10, 11},
                    {9, 10, 11, 12, 13, 14, 15},
                    {13, 14, 15, 16, 17, 18, 19},
                    {0, 1, 2, 3, 4, 5, 6},
                },
                new int[,] {
                    {1},
                    {5},
                    {9},
                    {13},
                },
                new int[,] {
                    {1, 2, 3, 4, 5, 6},
                    {5, 6, 7, 8, 9, 11},
                    {0, 2, 3, 4, 5, 6},
                },
                new int[,] {
                    {1, 2, 4},
                    {5, -6, 7},
                    {-1, 4, -5},
                    {1, 4, 5},
                },
                new int[,] {
                    {1, 2, 3, 4, 5, 6, 7},
                    {5, 6, 7, 8, -9, 10, 11},
                    {9, 10, -11, -12, -13, -14, -15},
                    {-13, -14, 15, 16, 17, 18, -19},
                },
                new int[,] {
                    {1, 3},
                    {5, -17},
                    {0, -3},
                },
                new int[,] {
                    {-9, -10, -11, -14, -15, 6},
                },
                new int[,] {
                    {1, 2, 3, 4, -5, -6, -7},
                    {5, 11, -17, 11, -10, 6, 5},
                    {-9, -10, -11, -14, -15, -16, 1},
                    {-9, -10, -11, -14, -15, -6, -2},
                    {-9, -10, -11, -14, -15, 6, 4},
                },
                new int[,] {
                    {1, 2, 3, 4, -5, -6, -7},
                    {5, 11, -17, 11, -10, 6, 5},
                    {-9, -10, -11, -14, -15, -16, 1},
                    {-9, -10, -11, -14, 15, -6, -2},
                    {-9, -10, -11, -14, -15, 6, 4},
                    {5, 11, -17, 11, -10, 6, -5},
                    {1, 1, -2, 3, -4, 0, 0},
                    {0, -2, -3, -4, -5, 0, 5},
                },
                new int[,] {
                    {1, 3, 4, -5},
                    {5, -17, 11, 7},
                    {-9, -11, -14, -15},
                    {-9, -11, -14, -6},
                    {0, -3, -4, -5},
                }
            };
            // Act
            for (int i = 0; i < answer.Length; i++)
            {
                _main.Task3(ref input[i], _main.FindLowerColIndex);
            }
            // Assert
            for (int i = 0; i < answer.Length; i++)
            {
                Assert.AreEqual(answer[i].GetLength(0), input[i].GetLength(0));
                for (int j = 0; j < answer[i].GetLength(0); j++)
                {
                    Assert.AreEqual(answer[i].GetLength(1), input[i].GetLength(1));
                    for (int k = 0; k < answer[i].GetLength(1); k++)
                    {
                        Assert.AreEqual(answer[i][j, k], input[i][j, k]);
                    }
                }
            }
        }
        [TestMethod]
        public void Test04()
        {
            // Arrange
            var input = _data.GetMatrixes();
            var answer = new int[][,] {
                new int[,] {
                    {1},
                    {5},
                    {9},
                    {13},
                    {0},
                },
                new int[,] {
                    {},
                    {},
                    {},
                    {},
                },
                new int[,] {
                    {1},
                    {5},
                    {0},
                },
                new int[,] {
                    {},
                    {},
                    {},
                    {},
                },
                new int[,] {
                    {},
                    {},
                    {},
                    {},
                },
                new int[,] {
                    {1},
                    {5},
                    {0},
                },
                new int[,] {
                    {},
                },
                new int[,] {
                    {},
                    {},
                    {},
                    {},
                    {},
                },
                new int[,] {
                    {1, -6, -7},
                    {5, 6, 5},
                    {-9, -16, 1},
                    {-9, -6, -2},
                    {-9, 6, 4},
                    {5, 6, -5},
                    {1, 0, 0},
                    {0, 0, 5},
                },
                new int[,] {
                    {1},
                    {5},
                    {-9},
                    {-9},
                    {0},
                }
            };
            // Act
            for (int i = 0; i < answer.Length; i++)
            {
                _main.Task4(ref input[i]);
            }
            // Assert
            for (int i = 0; i < answer.Length; i++)
            {
                Assert.AreEqual(answer[i].GetLength(0), input[i].GetLength(0));
                for (int j = 0; j < answer[i].GetLength(0); j++)
                {
                    Assert.AreEqual(answer[i].GetLength(1), input[i].GetLength(1));
                    for (int k = 0; k < answer[i].GetLength(1); k++)
                    {
                        Assert.AreEqual(answer[i][j, k], input[i][j, k]);
                    }
                }
            }
        }
        [TestMethod]
        public void Test05A()
        {
            // Arrange
            var input = _data.GetMatrixes();
            var answer = new int[][,] {
                new int[,] {
                    {1, 2, 3, 4, 5, 6, 7},
                    {5, 6, 7, 8, 9, 10, 11},
                    {9, 10, 11, 12, 13, 14, 15},
                    {0, 1, 2, 3, 4, 5, 6},
                },
                new int[,] {
                    {1},
                    {5},
                    {9},
                },
                new int[,] {
                    {1, 2, 3, 4, 5, 6},
                    {0, 2, 3, 4, 5, 6},
                },
                new int[,] {
                    {1, 2, 4, 6},
                    {-1, 4, -5, 6},
                    {1, 4, 5, 6},
                },
                new int[,] {
                    {1, 2, 3, 4, 5, 6, 7},
                    {5, 6, 7, 8, -9, 10, 11},
                    {9, 10, -11, -12, -13, -14, -15},
                },
                new int[,] {
                    {1, 2, 3},
                    {0, -2, -3},
                },
                new int[,] {
                },
                new int[,] {
                    {1, 2, 3, 4, -5, -6, -7},
                    {-9, -10, -11, -14, -15, -16, 1},
                    {-9, -10, -11, -14, -15, -6, -2},
                    {-9, -10, -11, -14, -15, 6, 4},
                },
                new int[,] {
                    {1, 2, 3, 4, -5, -6, -7},
                    {5, 11, -17, 11, -10, 6, 5},
                    {-9, -10, -11, -14, -15, -16, 1},
                    {-9, -10, -11, -14, -15, 6, 4},
                    {5, 11, -17, 11, -10, 6, -5},
                    {1, 1, -2, 3, -4, 0, 0},
                    {0, -2, -3, -4, -5, 0, 5},
                },
                new int[,] {
                    {1, 2, 3, 4, -5},
                    {-9, -10, -11, -14, -15},
                    {-9, -10, -11, -14, -6},
                    {0, -2, -3, -4, -5},
                }
            };
            // Act
            for (int i = 0; i < answer.Length; i++)
            {
                _main.Task5(ref input[i], _main.FindMax);
            }
            // Assert
            for (int i = 0; i < answer.Length; i++)
            {
                Assert.AreEqual(answer[i].GetLength(0), input[i].GetLength(0));
                for (int j = 0; j < answer[i].GetLength(0); j++)
                {
                    Assert.AreEqual(answer[i].GetLength(1), input[i].GetLength(1));
                    for (int k = 0; k < answer[i].GetLength(1); k++)
                    {
                        Assert.AreEqual(answer[i][j, k], input[i][j, k]);
                    }
                }
            }
        }
        [TestMethod]
        public void Test05B()
        {
            // Arrange
            var input = _data.GetMatrixes();
            var answer = new int[][,] {
                new int[,] {
                    {1, 2, 3, 4, 5, 6, 7},
                    {5, 6, 7, 8, 9, 10, 11},
                    {9, 10, 11, 12, 13, 14, 15},
                    {13, 14, 15, 16, 17, 18, 19},
                },
                new int[,] {
                    {5},
                    {9},
                    {13},
                },
                new int[,] {
                    {1, 2, 3, 4, 5, 6},
                    {5, 6, 7, 8, 9, 11},
                },
                new int[,] {
                    {1, 2, 4, 6},
                    {-1, 4, -5, 6},
                    {1, 4, 5, 6},
                },
                new int[,] {
                    {1, 2, 3, 4, 5, 6, 7},
                    {5, 6, 7, 8, -9, 10, 11},
                    {9, 10, -11, -12, -13, -14, -15},
                },
                new int[,] {
                    {1, 2, 3},
                    {0, -2, -3},
                },
                new int[,] {
                },
                new int[,] {
                    {1, 2, 3, 4, -5, -6, -7},
                    {-9, -10, -11, -14, -15, -16, 1},
                    {-9, -10, -11, -14, -15, -6, -2},
                    {-9, -10, -11, -14, -15, 6, 4},
                },
                new int[,] {
                    {1, 2, 3, 4, -5, -6, -7},
                    {-9, -10, -11, -14, -15, -16, 1},
                    {-9, -10, -11, -14, 15, -6, -2},
                    {-9, -10, -11, -14, -15, 6, 4},
                    {1, 1, -2, 3, -4, 0, 0},
                    {0, -2, -3, -4, -5, 0, 5},
                },
                new int[,] {
                    {1, 2, 3, 4, -5},
                    {-9, -10, -11, -14, -15},
                    {-9, -10, -11, -14, -6},
                    {0, -2, -3, -4, -5},
                }
            };
            // Act
            for (int i = 0; i < answer.Length; i++)
            {
                _main.Task5(ref input[i], _main.FindMin);
            }
            // Assert
            for (int i = 0; i < answer.Length; i++)
            {
                Assert.AreEqual(answer[i].GetLength(0), input[i].GetLength(0));
                for (int j = 0; j < answer[i].GetLength(0); j++)
                {
                    Assert.AreEqual(answer[i].GetLength(1), input[i].GetLength(1));
                    for (int k = 0; k < answer[i].GetLength(1); k++)
                    {
                        Assert.AreEqual(answer[i][j, k], input[i][j, k]);
                    }
                }
            }
        }
        [TestMethod]
        public void Test06A()
        {
            // Arrange
            var input = _data.GetMatrixes();
            var answer = new int[][,] {
                new int[,] {
                    {1, 2, 3, 4, 5, 6, 7},
                    {5, 6, 7, 8, 9, 10, 11},
                    {9, 10, 11, 12, 13, 14, 15},
                    {13, 14, 15, 16, 17, 18, 19},
                    {0, 1, 2, 3, 4, 5, 6},
                },
                new int[,] {
                    {1},
                    {5},
                    {9},
                    {13},
                },
                new int[,] {
                    {1, 2, 3, 4, 5, 6},
                    {5, 6, 7, 8, 9, 11},
                    {0, 2, 3, 4, 5, 6},
                },
                new int[,] {
                    {1, 2, 4, 6},
                    {5, -6, 7, 11},
                    {-1, 4, -5, 6},
                    {1, 4, 5, 6},
                },
                new int[,] {
                    {1, 2, 3, 4, 5, 6, 7},
                    {5, 6, 7, 8, -9, 10, 11},
                    {9, 10, -11, -12, -13, -14, -15},
                    {-19, -14, -13, 15, 16, 17, 18},
                },
                new int[,] {
                    {1, 2, 3},
                    {5, 11, -17},
                    {0, -2, -3},
                },
                new int[,] {
                    {-15, -14, -11, -10, -9, 6},
                },
                new int[,] {
                    {-7, -6, -5, 1, 2, 3, 4},
                    {5, 11, -17, 11, -10, 6, 5},
                    {-9, -10, -11, -14, -15, -16, 1},
                    {-15, -14, -11, -10, -9, -6, -2},
                    {-9, -10, -11, -14, -15, 6, 4},
                },
                new int[,] {
                    {-7, -6, -5, 1, 2, 3, 4},
                    {5, 11, -17, 11, -10, 6, 5},
                    {-9, -10, -11, -14, -15, -16, 1},
                    {-14, -11, -10, -9, -6, -2, 15},
                    {-9, -10, -11, -14, -15, 6, 4},
                    {5, 11, -17, 11, -10, 6, -5},
                    {-4, -2, 0, 0, 1, 1, 3},
                    {0, -2, -3, -4, -5, 0, 5},
                },
                new int[,] {
                    {-5, 1, 2, 3, 4},
                    {5, 11, -17, 11, 7},
                    {-9, -10, -11, -14, -15},
                    {-14, -11, -10, -9, -6},
                    {0, -2, -3, -4, -5},
                }
            };
            // Act
            for (int i = 0; i < answer.Length; i++)
            {
                _main.Task6(input[i], _main.SortRowAscending);
            }
            // Assert
            for (int i = 0; i < answer.Length; i++)
            {
                Assert.AreEqual(answer[i].GetLength(0), input[i].GetLength(0));
                for (int j = 0; j < answer[i].GetLength(0); j++)
                {
                    Assert.AreEqual(answer[i].GetLength(1), input[i].GetLength(1));
                    for (int k = 0; k < answer[i].GetLength(1); k++)
                    {
                        Assert.AreEqual(answer[i][j, k], input[i][j, k]);
                    }
                }
            }
        }

        [TestMethod]
        public void Test06B()
        {
            // Arrange
            var input = _data.GetMatrixes();
            var answer = new int[][,] {
                new int[,] {
                    {7, 6, 5, 4, 3, 2, 1},
                    {5, 6, 7, 8, 9, 10, 11},
                    {9, 10, 11, 12, 13, 14, 15},
                    {19, 18, 17, 16, 15, 14, 13},
                    {0, 1, 2, 3, 4, 5, 6},
                },
                new int[,] {
                    {1},
                    {5},
                    {9},
                    {13},
                },
                new int[,] {
                    {6, 5, 4, 3, 2, 1},
                    {5, 6, 7, 8, 9, 11},
                    {0, 2, 3, 4, 5, 6},
                },
                new int[,] {
                    {6, 4, 2, 1},
                    {5, -6, 7, 11},
                    {-1, 4, -5, 6},
                    {6, 5, 4, 1},
                },
                new int[,] {
                    {7, 6, 5, 4, 3, 2, 1},
                    {5, 6, 7, 8, -9, 10, 11},
                    {9, 10, -11, -12, -13, -14, -15},
                    {18, 17, 16, 15, -13, -14, -19},
                },
                new int[,] {
                    {3, 2, 1},
                    {5, 11, -17},
                    {0, -2, -3},
                },
                new int[,] {
                    {6, -9, -10, -11, -14, -15},
                },
                new int[,] {
                    {4, 3, 2, 1, -5, -6, -7},
                    {5, 11, -17, 11, -10, 6, 5},
                    {-9, -10, -11, -14, -15, -16, 1},
                    {-2, -6, -9, -10, -11, -14, -15},
                    {-9, -10, -11, -14, -15, 6, 4},
                },
                new int[,] {
                    {4, 3, 2, 1, -5, -6, -7},
                    {5, 11, -17, 11, -10, 6, 5},
                    {-9, -10, -11, -14, -15, -16, 1},
                    {15, -2, -6, -9, -10, -11, -14},
                    {-9, -10, -11, -14, -15, 6, 4},
                    {5, 11, -17, 11, -10, 6, -5},
                    {3, 1, 1, 0, 0, -2, -4},
                    {0, -2, -3, -4, -5, 0, 5},
                },
                new int[,] {
                    {4, 3, 2, 1, -5},
                    {5, 11, -17, 11, 7},
                    {-9, -10, -11, -14, -15},
                    {-6, -9, -10, -11, -14},
                    {0, -2, -3, -4, -5},
                }
            };
            // Act
            for (int i = 0; i < answer.Length; i++)
            {
                _main.Task6(input[i], _main.SortRowDescending);
            }
            // Assert
            for (int i = 0; i < answer.Length; i++)
            {
                Assert.AreEqual(answer[i].GetLength(0), input[i].GetLength(0));
                for (int j = 0; j < answer[i].GetLength(0); j++)
                {
                    Assert.AreEqual(answer[i].GetLength(1), input[i].GetLength(1));
                    for (int k = 0; k < answer[i].GetLength(1); k++)
                    {
                        Assert.AreEqual(answer[i][j, k], input[i][j, k]);
                    }
                }
            }
        }
        [TestMethod]
        public void Test07A()
        {
            // Arrange
            var input = _data.GetMatrixes();
            var answer = new int[][,] {
                new int[,] {
                    {1, 2, 3, 4, 5, 6, 0},
                    {5, 6, 7, 8, 9, 10, 0},
                    {9, 10, 11, 12, 13, 14, 0},
                    {13, 14, 15, 16, 17, 18, 0},
                    {0, 1, 2, 3, 4, 5, 0},
                },
                new int[,] {
                    {0},
                    {0},
                    {0},
                    {0},
                },
                new int[,] {
                    {1, 2, 3, 4, 5, 0},
                    {5, 6, 7, 8, 9, 0},
                    {0, 2, 3, 4, 5, 0},
                },
                new int[,] {
                    {1, 2, 4, 0},
                    {5, -6, 7, 0},
                    {-1, 4, -5, 0},
                    {1, 4, 5, 0},
                },
                new int[,] {
                    {1, 2, 3, 4, 5, 6, 0},
                    {5, 6, 7, 8, -9, 10, 0},
                    {9, 0, -11, -12, -13, -14, -15},
                    {-13, -14, 15, 16, 17, 0, -19},
                },
                new int[,] {
                    {1, 2, 0},
                    {5, 0, -17},
                    {0, -2, -3},
                },
                new int[,] {
                    {-9, -10, -11, -14, -15, 0},
                },
                new int[,] {
                    {1, 2, 3, 0, -5, -6, -7},
                    {5, 0, -17, 0, -10, 6, 5},
                    {-9, -10, -11, -14, -15, -16, 0},
                    {-9, -10, -11, -14, -15, -6, 0},
                    {-9, -10, -11, -14, -15, 0, 4},
                },
                new int[,] {
                    {1, 2, 3, 0, -5, -6, -7},
                    {5, 0, -17, 0, -10, 6, 5},
                    {-9, -10, -11, -14, -15, -16, 0},
                    {-9, -10, -11, -14, 0, -6, -2},
                    {-9, -10, -11, -14, -15, 0, 4},
                    {5, 0, -17, 0, -10, 6, -5},
                    {1, 1, -2, 0, -4, 0, 0},
                    {0, -2, -3, -4, -5, 0, 0},
                },
                new int[,] {
                    {1, 2, 3, 0, -5},
                    {5, 0, -17, 0, 7},
                    {0, -10, -11, -14, -15},
                    {-9, -10, -11, -14, 0},
                    {0, -2, -3, -4, -5},
                }
            };
            // Act
            for (int i = 0; i < answer.Length; i++)
            {
                _main.Task7(input[i], _main.ReplaceByZero);
            }
            // Assert
            for (int i = 0; i < answer.Length; i++)
            {
                Assert.AreEqual(answer[i].GetLength(0), input[i].GetLength(0));
                for (int j = 0; j < answer[i].GetLength(0); j++)
                {
                    Assert.AreEqual(answer[i].GetLength(1), input[i].GetLength(1));
                    for (int k = 0; k < answer[i].GetLength(1); k++)
                    {
                        Assert.AreEqual(answer[i][j, k], input[i][j, k]);
                    }
                }
            }
        }
        [TestMethod]
        public void Test07B()
        {
            // Arrange
            var input = _data.GetMatrixes();
            var answer = new int[][,] {
                new int[,] {
                    {1, 2, 3, 4, 5, 6, 49},
                    {5, 6, 7, 8, 9, 10, 77},
                    {9, 10, 11, 12, 13, 14, 105},
                    {13, 14, 15, 16, 17, 18, 133},
                    {0, 1, 2, 3, 4, 5, 42},
                },
                new int[,] {
                    {1},
                    {5},
                    {9},
                    {13},
                },
                new int[,] {
                    {1, 2, 3, 4, 5, 36},
                    {5, 6, 7, 8, 9, 66},
                    {0, 2, 3, 4, 5, 36},
                },
                new int[,] {
                    {1, 2, 4, 24},
                    {5, -6, 7, 44},
                    {-1, 4, -5, 24},
                    {1, 4, 5, 24},
                },
                new int[,] {
                    {1, 2, 3, 4, 5, 6, 49},
                    {5, 6, 7, 8, -9, 10, 77},
                    {9, 20, -11, -12, -13, -14, -15},
                    {-13, -14, 15, 16, 17, 108, -19},
                },
                new int[,] {
                    {1, 2, 9},
                    {5, 22, -17},
                    {0, -2, -3},
                },
                new int[,] {
                    {-9, -10, -11, -14, -15, 36},
                },
                new int[,] {
                    {1, 2, 3, 16, -5, -6, -7},
                    {5, 22, -17, 44, -10, 6, 5},
                    {-9, -10, -11, -14, -15, -16, 7},
                    {-9, -10, -11, -14, -15, -6, -14},
                    {-9, -10, -11, -14, -15, 36, 4},
                },
                new int[,] {
                    {1, 2, 3, 16, -5, -6, -7},
                    {5, 22, -17, 44, -10, 6, 5},
                    {-9, -10, -11, -14, -15, -16, 7},
                    {-9, -10, -11, -14, 75, -6, -2},
                    {-9, -10, -11, -14, -15, 36, 4},
                    {5, 22, -17, 44, -10, 6, -5},
                    {1, 1, -2, 12, -4, 0, 0},
                    {0, -2, -3, -4, -5, 0, 35},
                },
                new int[,] {
                    {1, 2, 3, 16, -5},
                    {5, 22, -17, 44, 7},
                    {-9, -10, -11, -14, -15},
                    {-9, -10, -11, -14, -30},
                    {0, -2, -3, -4, -5},
                }
            };
            // Act
            for (int i = 0; i < answer.Length; i++)
            {
                _main.Task7(input[i], _main.MultiplyByColumn);
            }
            // Assert
            for (int i = 0; i < answer.Length; i++)
            {
                Assert.AreEqual(answer[i].GetLength(0), input[i].GetLength(0));
                for (int j = 0; j < answer[i].GetLength(0); j++)
                {
                    Assert.AreEqual(answer[i].GetLength(1), input[i].GetLength(1));
                    for (int k = 0; k < answer[i].GetLength(1); k++)
                    {
                        Assert.AreEqual(answer[i][j, k], input[i][j, k]);
                    }
                }
            }
        }
        [TestMethod]
        public void Test08A()
        {
            // Arrange
            var inputA = 0.1;
            var inputB = 1;
            var inputH = 0.1;
            var answer = new double[,] {
                    {2.6912681270017194, 2.691268139166703},
                    {2.6122212502032887, 2.6122204929844544},
                    {2.4868568953361008, 2.486856868603152},
                    {2.323884253225943, 2.3238842457941966},
                    {2.1339300133437806, 2.133930111437405},
                    {1.9283342127279457, 1.9283342378052784},
                    {1.7179999583568928, 1.7179999609519054},
                    {1.5124670271786733, 1.5124670047163078},
                    {1.3193029843812667, 1.3193027107322826},
                    {1.1438356417724056, 1.1438356437916406},
                };
            var test = new double[answer.GetLength(0), answer.GetLength(1)];
            // Act
            test = _main.Task8(inputA, inputB, inputH, _main.SumA, _main.YA);
            // Assert
            Assert.AreEqual(answer.GetLength(0), test.GetLength(0));
            Assert.AreEqual(answer.GetLength(1), test.GetLength(1));
            for (int i = 0; i < answer.GetLength(0); i++)
            {
                for (int j = 0; j < answer.GetLength(1); j++)
                {
                    Assert.AreEqual(answer[i, j], test[i, j], E);
                }
            }
        }
        [TestMethod]
        public void Test08B()
        {
            // Arrange
            var inputA = Math.PI / 5;
            var inputB = Math.PI;
            var inputH = Math.PI / 25;
            var answer = new double[,] {
                    {-7.303507257791954, -7.303507256806125},
                    {-7.260083616438886, -7.260080997441332},
                    {-7.20875508832593, -7.208759054555667},
                    {-7.149545913780194, -7.149541428149131},
                    {-7.082422850960618, -7.082428118221723},
                    {-7.007419743934666, -7.0074191247734445},
                    {-6.924507777616338, -6.9245144478042935},
                    {-6.833707144188682, -6.833714087314272},
                    {-6.735009114927305, -6.735018043303378},
                    {-6.628434910121376, -6.6284263157716135},
                    {-6.51393629032441, -6.513938904718977},
                    {-6.391544374889424, -6.391555810145468},
                    {-6.261289902367042, -6.261277032051089},
                    {-6.123087221829428, -6.1231025704358375},
                    {-5.977012101846995, -5.977032425299715},
                    {-5.823071802317092, -5.823066596642721},
                    {-5.661174273430562, -5.661205084464855},
                    {-5.49140693658761, -5.491447888766118},
                    {-5.313734462977021, -5.313795009546508},
                    {-5.128131509726421, -5.128246446806029},
                    {-4.935801700711345, -4.934802200544677},
                };
            var test = new double[answer.GetLength(0), answer.GetLength(1)];
            // Act
            test = _main.Task8(inputA, inputB, inputH, _main.SumB, _main.YB);
            // Assert
            Assert.AreEqual(answer.GetLength(0), test.GetLength(0));
            Assert.AreEqual(answer.GetLength(1), test.GetLength(1));
            for (int i = 0; i < answer.GetLength(0); i++)
            {
                for (int j = 0; j < answer.GetLength(1); j++)
                {
                    Assert.AreEqual(answer[i, j], test[i, j], E);
                }
            }
        }
        [TestMethod]
        public void Test09A()
        {
            // Arrange
            var input = _data.GetMatrixes();
            var answer = new int[] { 0, 0, 0, 360, 0, 433, 0, 0, 0, 1434 };
            var test = new int[answer.Length];
            // Act
            for (int i = 0; i < answer.Length; i++)
            {
                test[i] = _main.Task9(input[i], _main.GetUpperTriangle);
            }
            // Assert
            Assert.AreEqual(answer.Length, test.Length);
            for (int i = 0; i < answer.Length; i++)
            {
                Assert.AreEqual(answer[i], test[i]);
            }
        }
        [TestMethod]
        public void Test09B()
        {
            // Arrange
            var input = _data.GetMatrixes();
            var answer = new int[] { 0, 0, 0, 182, 0, 160, 0, 0, 0, 1001 };
            var test = new int[answer.Length];
            // Act
            for (int i = 0; i < answer.Length; i++)
            {
                test[i] = _main.Task9(input[i], _main.GetLowerTriangle);
            }
            // Assert
            Assert.AreEqual(answer.Length, test.Length);
            for (int i = 0; i < answer.Length; i++)
            {
                Assert.AreEqual(answer[i], test[i]);
            }
        }
        [TestMethod]
        public void Test10A()
        {
            // Arrange
            var input = _data.GetArrayArrays();
            var answer = new bool[5] { false, true, true, true, false };
            var test = new bool[answer.Length];
            // Act
            for (int i = 0; i < answer.Length; i++)
            {
                test[i] = _main.Task10(input[i], _main.CheckTransformAbility);
            }
            // Assert
            Assert.AreEqual(answer.Length, test.Length);
            for (int i = 0; i < answer.Length; i++)
            {
                Assert.AreEqual(answer[i], test[i]);
            }
        }
        [TestMethod]
        public void Test10B()
        {
            // Arrange
            var input = _data.GetArrayArrays();
            var answer = new bool[5] { false, false, true, false, false };
            var test = new bool[answer.Length];
            // Act
            for (int i = 0; i < answer.Length; i++)
            {
                test[i] = _main.Task10(input[i], _main.CheckSumOrder);
            }
            // Assert
            Assert.AreEqual(answer.Length, test.Length);
            for (int i = 0; i < answer.Length; i++)
            {
                Assert.AreEqual(answer[i], test[i]);
            }
        }
        [TestMethod]
        public void Test10C()
        {
            // Arrange
            var input = _data.GetArrayArrays();
            var answer = new bool[5] { true, true, false, true, false };
            var test = new bool[answer.Length];
            // Act
            for (int i = 0; i < answer.Length; i++)
            {
                test[i] = _main.Task10(input[i], _main.CheckArraysOrder);
            }
            // Assert
            Assert.AreEqual(answer.Length, test.Length);
            for (int i = 0; i < answer.Length; i++)
            {
                Assert.AreEqual(answer[i], test[i]);
            }
        }
        [TestMethod]
        public void Test_FindDiagonalMaxIndex()
        {
            var matrix = new int[,] { { 1, 2 }, { 3, 4 } };
            int result = _main.FindDiagonalMaxIndex(matrix);
            Assert.AreEqual(1, result);
        }

        [TestMethod]
        public void Test_RemoveRow()
        {
            var matrix = new int[,] { { 1, 2 }, { 3, 4 }, { 5, 6 } };
            _main.RemoveRow(ref matrix, 1);
            var expected = new int[,] { { 1, 2 }, { 5, 6 } };
            CollectionAssert.AreEqual(expected.Cast<int>().ToArray(), matrix.Cast<int>().ToArray());
        }

        [TestMethod]
        public void Test_GetAverageExceptEdges()
        {
            var matrix = new int[,] { { 1, 2, 3 }, { 4, 5, 6 } };
            double avg = _main.GetAverageExceptEdges(matrix);
            // Убираем 1 и 6 → остаются [2,3,4,5], среднее = 3.5
            Assert.AreEqual(3.5, avg, E);
        }

        [TestMethod]
        public void Test_FindMaxMin()
        {
            var matrix = new int[,] { { -2, 10 }, { -5, 0 } };
            Assert.AreEqual(10, _main.FindMax(matrix, out int i, out int j));
            Assert.AreEqual(0, i);
            Assert.AreEqual(1, j);
            Assert.AreEqual(-5, _main.FindMin(matrix, out i, out j));
            Assert.AreEqual(1, i);
            Assert.AreEqual(0, j);
        }

        [TestMethod]
        public void Test_RemoveColumn()
        {
            var matrix = new int[,] { { 1, 2, 3 }, { 4, 5, 6 } };
            _main.RemoveColumn(ref matrix, 1);
            var expected = new int[,] { { 1, 3 }, { 4, 6 } };
            CollectionAssert.AreEqual(expected.Cast<int>().ToArray(), matrix.Cast<int>().ToArray());
        }

        [TestMethod]
        public void Test_CheckZerosInColumn()
        {
            var matrix = new int[,] { { 1, 0 }, { 3, 4 } };
            Assert.IsTrue(_main.CheckZerosInColumn(matrix, 1));
            Assert.IsFalse(_main.CheckZerosInColumn(matrix, 0));
        }

        [TestMethod]
        public void Test_FindMaxMinWithOut()
        {
            var matrix = new int[,] { { 1, 2 }, { 3, 4 } };
            int row, col;
            Assert.AreEqual(4, _main.FindMax(matrix, out row, out col));
            Assert.AreEqual(1, _main.FindMin(matrix, out row, out col));
        }

        [TestMethod]
        public void Test_SortRowAscendingDescending()
        {
            var matrix = new int[,] { { 3, 1, 2 },
            { 8, 2, 54 } };
            _main.SortRowAscending(matrix, 0);
            CollectionAssert.AreEqual(new int[] { 1, 2, 3, 8, 2, 54 }, matrix.Cast<int>().ToArray());

            matrix = new int[,] { { 3, 1, 2 },
            { 8, 2, 54 } };
            _main.SortRowDescending(matrix, 0);
            CollectionAssert.AreEqual(new int[] { 3, 2, 1, 8, 2, 54 }, matrix.Cast<int>().ToArray());
        }

        [TestMethod]
        public void Test_FindMaxInRowAndReplace()
        {
            var matrix = new int[,] { { 1, 3, 3 } };
            int max = _main.FindMaxInRow(matrix, 0);
            Assert.AreEqual(3, max);

            _main.ReplaceByZero(matrix, 0, max);
            CollectionAssert.AreEqual(new int[] { 1, 0, 0 }, matrix.Cast<int>().ToArray());
        }

        [TestMethod]
        public void Test_MultiplyByColumn()
        {
            var matrix = new int[,] { { 2, 5, 5 } };
            _main.MultiplyByColumn(matrix, 0, 5);
            // 5*2 и 5*3 → 10, 15
            CollectionAssert.AreEqual(new int[] { 2, 10, 15 }, matrix.Cast<int>().ToArray());
        }

        [TestMethod]
        public void Test_Sum_Simple()
        {
            var arr = new int[] { 1, 2, 3 };
            int result = _main.Sum(arr);
            // 1²+2²+3²=14
            Assert.AreEqual(14, result);
        }

        [TestMethod]
        public void Test_GetUpperAndLowerTriangle()
        {
            var matrix = new int[,] { { 1, 2 }, { 3, 4 } };
            var upper = _main.GetUpperTriangle(matrix);
            var lower = _main.GetLowerTriangle(matrix);
            CollectionAssert.AreEqual(new int[] { 1, 2, 4 }, upper);
            CollectionAssert.AreEqual(new int[] { 1, 3, 4 }, lower);
        }

        [TestMethod]
        public void Test_CheckSumOrder()
        {
            var array = new int[][] { new[] { 1, 2 }, new[] { 3, 4 }, new[] { 5, 6 } };
            Assert.IsTrue(_main.CheckSumOrder(array)); // возрастающая сумма
        }

        [TestMethod]
        public void Test_CheckArraysOrder()
        {
            var array = new int[][] { new[] { 1, 2, 3 }, new[] { 3, 1 } };
            Assert.IsTrue(_main.CheckArraysOrder(array)); // один из массивов упорядочен
        }

        [TestMethod]
        public void Test_CheckTransformAbility()
        {
            var array = new int[][] { new[] { 1, 2, 3 }, new[] { 4, 5 }, new[] { 6 } };
            Assert.IsTrue(_main.CheckTransformAbility(array));
        }
    }
}
