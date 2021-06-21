using System.Collections;
using System.Collections.Generic;
using UnityEngine;

public class MySparseMatrix
{
    public List<(int col, float v)>[] data;
    int n;
    public MySparseMatrix(int _n)
    {
        n = _n;
        data = new List<(int col, float v)>[n];
        for (int i = 0; i < n; i++)
        {
            data[i] = new List<(int col, float v)>();
        }
    }
    public void Insert(int row, int col, float v)
    {
        data[row].Add((col, v));
    }
    public void Modify(int row, int col, float v)
    {
        for (int i = 0; i < data[row].Count; i++)
        {
            if (data[row][i].col == col)
            {
                data[row][i] = (col, v);
                return;
            }
        }
        Insert(row, col, v);
    }
    public void Add(int row, int col, float v)
    {
        for (int i = 0; i < data[row].Count; i++)
        {
            if (data[row][i].col == col)
            {
                data[row][i] = (col, data[row][i].v + v);
                return;
            }
        }
        Insert(row, col, v);
    }
    float read(int row, int col)
    {
        for (int i = 0; i < data[row].Count; i++)
        {
            if (data[row][i].col == col)
            {
                return data[row][i].v;
            }
        }
        return 0;
    }
    public void Multiply(ref float[] x, ref float[] y) {
        for (int i = 0; i < n; i++)
        {
            y[i] = 0;
            for (int j = 0; j < data[i].Count; j++)
            {
                y[i] += x[data[i][j].col] * data[i][j].v;
            }
        }
    }
    public void MultiplyByScalar(float x)
    {
        for (int i = 0; i < n; i++)
        {
            for (int j = 0; j < data[i].Count; j++)
            {
                data[i][j] = (data[i][j].col, data[i][j].v * x);
            }
        }
    }
}
