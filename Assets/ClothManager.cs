using System.Collections;
using System.Collections.Generic;
using UnityEngine;
using UnityEditor;

public class ClothManager : MonoBehaviour
{
    // Start is called before the first frame update
    static public Mesh terrainMesh;
    const float Size = 10f;
    [Range(3, 100)] public int N = 10;
    public float totMass = 100;
    Vector3[] f;
    Vector3[] v;
    float[] v0;
    float[] f0;
    float[] B;
    float[] X;
    float dT = 0.03f;
    Matrix4x4 renderMatrix = new Matrix4x4();
    List<Mesh> grassMeshes = new List<Mesh>();
    public float kStructuralC = 10f;
    public float kBendC = 1f;
    public float kDamp = 5f;
    public float fWind = 1f;
    public float tWind = 1f;
    float tnow = 0;
    MySparseMatrix mat;
    void Start()
    {
        f = new Vector3[N * N + (N - 1) * (N - 1)];
        v = new Vector3[N * N + (N - 1) * (N - 1)];
        f0 = new float[f.Length * 3];
        v0 = new float[v.Length * 3];
        B = new float[v.Length * 3];
        X = new float[v.Length * 3];
        BuildClothMesh();
        Application.targetFrameRate = 60;
    }
    void CalcNormal(ref Vector3[] vertices, ref List<Vector3> normals)
    {
        for (int i = 0; i < N; i++)
        {
            for (int j = 0; j < N; j++)
            {
                Vector3 dx;
                Vector3 dz;
                if (i != N - 1 && i != 0)
                {
                    dx = vertices[i * N + j + N] - vertices[i * N + j - N];
                }
                else if (i == 0)
                {
                    dx = (vertices[i * N + j + N] - vertices[i * N + j]) * 2;
                }
                else
                {
                    dx = (vertices[i * N + j] - vertices[i * N + j - N]) * 2;
                }
                if (j != N - 1 && j != 0)
                {
                    dz = (vertices[i * N + j + 1] - vertices[i * N + j - 1]);
                }
                else if (j == 0)
                {
                    dz = (vertices[i * N + j + 1] - vertices[i * N + j]) * 2;
                }
                else
                {
                    dz = (vertices[i * N + j] - vertices[i * N + j - 1]) * 2;
                }
                normals.Add(Vector3.Cross(dz, dx).normalized);

            }
        }
        for (int i = 0; i < N - 1; i++)
        {
            for (int j = 0; j < N - 1; j++)
            {
                Vector3 dx;
                Vector3 dz;
                if (i != N - 2 && i != 0)
                {
                    dx = vertices[N * N + i * (N - 1) + j + N - 1] - vertices[N * N + i * (N - 1) + j - (N - 1)];
                }
                else if (i == 0)
                {
                    dx = (vertices[N * N + i * (N - 1) + j + N - 1] - vertices[N * N + i * (N - 1) + j]) * 2;
                }
                else
                {
                    dx = (vertices[N * N + i * (N - 1) + j] - vertices[N * N + i * (N - 1) + j - (N - 1)]) * 2;
                }
                if (j != N - 2 && j != 0)
                {
                    dz = (vertices[N * N + i * (N - 1) + j + 1] - vertices[N * N + i * (N - 1) + j - 1]);
                }
                else if (j == 0)
                {
                    dz = (vertices[N * N + i * (N - 1) + j + 1] - vertices[N * N + i * (N - 1) + j]) * 2;
                }
                else
                {
                    dz = (vertices[N * N + i * (N - 1) + j] - vertices[N * N + i * (N - 1) + j - 1]) * 2;
                }
                normals.Add(Vector3.Cross(dz, dx).normalized);

            }
        }
    }

    public void BuildClothMesh()
    {
        MeshFilter meshFilter = gameObject.GetComponent<MeshFilter>();

        terrainMesh = new Mesh();
        List<Vector3> vertices = new List<Vector3>();
        List<Vector2> uvs = new List<Vector2>();
        List<Vector3> normals = new List<Vector3>();
        List<int> triangles = new List<int>();
        for (int i = 0; i < N; i++)
        {
            for (int j = 0; j < N; j++)
            {
                Vector2 uv = new Vector2(1f * i / (N - 1), 1f * j / (N - 1));
                uvs.Add(uv);
                Vector3 p = new Vector3(uv.x - .5f, 0, uv.y - .5f) * 2f;
                Vector3 v = new Vector3(uv.x * Size,
                                        0f,
                                        uv.y * Size);
                /*
                if (p.magnitude > 1)
                {
                    v = p.normalized / 2 + new Vector3(.5f, 0f, .5f);
                    v.y = 1 - p.magnitude;
                    v *= Size;
                    normals.Add(p.normalized);
                } else {
                    p.y = Mathf.Sqrt(1 - Mathf.Pow(p.magnitude, 2));
                    v.y = p.y * Size / 2f;
                    normals.Add(p);
                }
                */
                vertices.Add(v);
            }
        }
        for (int i = 0; i < N - 1; i++)
        {
            for (int j = 0; j < N - 1; j++)
            {
                Vector2 uv = new Vector2((i + .5f) / (N - 1), (j + .5f) / (N - 1));
                uvs.Add(uv);
                Vector3 v = new Vector3(uv.x * Size,
                                        0f,
                                        uv.y * Size);
                vertices.Add(v);
            }
        }
        for (int i = 0; i < N - 1; i++)
        {
            for (int j = 0; j < N - 1; j++)
            {
                triangles.Add(N * N + i * (N - 1) + j);
                triangles.Add(i * N + j + 1);
                triangles.Add(i * N + j + N + 1);
                triangles.Add(N * N + i * (N - 1) + j);
                triangles.Add(i * N + j + N + 1);
                triangles.Add(i * N + j + N);
                triangles.Add(N * N + i * (N - 1) + j);
                triangles.Add(i * N + j + N);
                triangles.Add(i * N + j);
                triangles.Add(N * N + i * (N - 1) + j);
                triangles.Add(i * N + j);
                triangles.Add(i * N + j + 1);
            }
        }

        terrainMesh.SetVertices(vertices);
        terrainMesh.SetUVs(0, uvs);
        terrainMesh.SetTriangles(triangles, 0);
        terrainMesh.SetNormals(normals);
        meshFilter.mesh = terrainMesh;
        for (int i = 0; i < v.Length; i++)
        {
            v[i] = new Vector3(0f, 0f, 0f);
        }
        for (int i = 0; i < v0.Length; i++)
        {
            v0[i] = 0f;
        }
    }


    public void ForwardEuler(ref Vector3[] vertices)
    {

        float lRestStructural = Size / (N - 1);
        float lRestShear = lRestStructural / Mathf.Sqrt(2);
        float lRestBend = lRestStructural * 2;

        float kStructural = kStructuralC / (1f / (N - 1));
        float kShear = kStructuralC / (1f / (N - 1) * Mathf.Sqrt(2f) / 2);
        float kBend = kBendC / (1f / (N - 1) * 2);
        for (int i = 0; i < N; i++)
        {
            for (int j = 0; j < N; j++)
            {
                f[i * N + j] = new Vector3(0f, 0f, 0f);
                if (i - 2 >= 0)
                {
                    Vector3 d = vertices[(i - 2) * N + j] - vertices[i * N + j];
                    f[i * N + j] += kBend * d.normalized * (d.magnitude - lRestBend);
                }
                if (j - 2 >= 0)
                {
                    Vector3 d = vertices[i * N + j - 2] - vertices[i * N + j];
                    f[i * N + j] += kBend * d.normalized * (d.magnitude - lRestBend);
                }
                if (i + 2 < N)
                {
                    Vector3 d = vertices[(i + 2) * N + j] - vertices[i * N + j];
                    f[i * N + j] += kBend * d.normalized * (d.magnitude - lRestBend);
                }
                if (j + 2 < N)
                {
                    Vector3 d = vertices[i * N + j + 2] - vertices[i * N + j];
                    f[i * N + j] += kBend * d.normalized * (d.magnitude - lRestBend);
                }

                if (i - 1 >= 0)
                {
                    Vector3 d = vertices[(i - 1) * N + j] - vertices[i * N + j];
                    f[i * N + j] += kStructural * d.normalized * (d.magnitude - lRestStructural);
                }
                if (j - 1 >= 0)
                {
                    Vector3 d = vertices[i * N + j - 1] - vertices[i * N + j];
                    f[i * N + j] += kStructural * d.normalized * (d.magnitude - lRestStructural);
                }
                if (i + 1 < N)
                {
                    Vector3 d = vertices[(i + 1) * N + j] - vertices[i * N + j];
                    f[i * N + j] += kStructural * d.normalized * (d.magnitude - lRestStructural);
                }
                if (j + 1 < N)
                {
                    Vector3 d = vertices[i * N + j + 1] - vertices[i * N + j];
                    f[i * N + j] += kStructural * d.normalized * (d.magnitude - lRestStructural);
                }

                if (i - 1 >= 0 && j - 1 >= 0)
                {
                    Vector3 d = vertices[N * N + (i - 1) * (N - 1) + j - 1] - vertices[i * N + j];
                    f[i * N + j] += kShear * d.normalized * (d.magnitude - lRestShear);
                }
                if (i - 1 >= 0 && j < N - 1)
                {
                    Vector3 d = vertices[N * N + (i - 1) * (N - 1) + j] - vertices[i * N + j];
                    f[i * N + j] += kShear * d.normalized * (d.magnitude - lRestShear);
                }
                if (i < N - 1 && j - 1 >= 0)
                {
                    Vector3 d = vertices[N * N + i * (N - 1) + j - 1] - vertices[i * N + j];
                    f[i * N + j] += kShear * d.normalized * (d.magnitude - lRestShear);
                }
                if (i < N - 1 && j < N - 1)
                {
                    Vector3 d = vertices[N * N + i * (N - 1) + j] - vertices[i * N + j];
                    f[i * N + j] += kShear * d.normalized * (d.magnitude - lRestShear);
                }
            }
        }
        for (int i = 0; i < N - 1; i++)
        {
            for (int j = 0; j < N - 1; j++)
            {
                f[N * N + i * (N - 1) + j] = new Vector3(0f, 0f, 0f);
                if (i - 2 >= 0)
                {
                    Vector3 d = vertices[N * N + (i - 2) * (N - 1) + j] - vertices[N * N + i * (N - 1) + j];
                    f[N * N + i * (N - 1) + j] += kBend * d.normalized * (d.magnitude - lRestBend);
                }
                if (j - 2 >= 0)
                {
                    Vector3 d = vertices[N * N + i * (N - 1) + j - 2] - vertices[N * N + i * (N - 1) + j];
                    f[N * N + i * (N - 1) + j] += kBend * d.normalized * (d.magnitude - lRestBend);
                }
                if (i + 2 < N - 1)
                {
                    Vector3 d = vertices[N * N + (i + 2) * (N - 1) + j] - vertices[N * N + i * (N - 1) + j];
                    f[N * N + i * (N - 1) + j] += kBend * d.normalized * (d.magnitude - lRestBend);
                }
                if (j + 2 < N - 1)
                {
                    Vector3 d = vertices[N * N + i * (N - 1) + j + 2] - vertices[N * N + i * (N - 1) + j];
                    f[N * N + i * (N - 1) + j] += kBend * d.normalized * (d.magnitude - lRestBend);
                }
                if (i - 1 >= 0)
                {
                    Vector3 d = vertices[N * N + (i - 1) * (N - 1) + j] - vertices[N * N + i * (N - 1) + j];
                    f[N * N + i * (N - 1) + j] += kStructural * d.normalized * (d.magnitude - lRestStructural);
                }
                if (j - 1 >= 0)
                {
                    Vector3 d = vertices[N * N + i * (N - 1) + j - 1] - vertices[N * N + i * (N - 1) + j];
                    f[N * N + i * (N - 1) + j] += kStructural * d.normalized * (d.magnitude - lRestStructural);
                }
                if (i + 1 < N - 1)
                {
                    Vector3 d = vertices[N * N + (i + 1) * (N - 1) + j] - vertices[N * N + i * (N - 1) + j];
                    f[N * N + i * (N - 1) + j] += kStructural * d.normalized * (d.magnitude - lRestStructural);
                }
                if (j + 1 < N - 1)
                {
                    Vector3 d = vertices[N * N + i * (N - 1) + j + 1] - vertices[N * N + i * (N - 1) + j];
                    f[N * N + i * (N - 1) + j] += kStructural * d.normalized * (d.magnitude - lRestStructural);
                }

                {
                    Vector3 d = vertices[i * N + j] - vertices[N * N + i * (N - 1) + j];
                    f[N * N + i * (N - 1) + j] += kShear * d.normalized * (d.magnitude - lRestShear);
                }
                {
                    Vector3 d = vertices[i * N + j + 1] - vertices[N * N + i * (N - 1) + j];
                    f[N * N + i * (N - 1) + j] += kShear * d.normalized * (d.magnitude - lRestShear);
                }
                {
                    Vector3 d = vertices[(i + 1) * N + j] - vertices[N * N + i * (N - 1) + j];
                    f[N * N + i * (N - 1) + j] += kShear * d.normalized * (d.magnitude - lRestShear);
                }
                {
                    Vector3 d = vertices[(i + 1) * N + j + 1] - vertices[N * N + i * (N - 1) + j];
                    f[N * N + i * (N - 1) + j] += kShear * d.normalized * (d.magnitude - lRestShear);
                }
            }
        }
        float mass = totMass / v.Length;
        for (int i = 0; i < f.Length; i++) {
            f[i] += new Vector3(0f, -1f, 0f) * 9.8f * mass;
        }
        for (int i = 0; i < f.Length; i++)
        {
            f[i] += -v[i] * kDamp;
        }

        tnow += dT;
        for (int i = 0; i < f.Length; i++)
        {
            f[i] += new Vector3(0f, 0f, -1f) * Mathf.Pow(Mathf.Sin(tnow * tWind), 2f) * fWind;
        }


        for (int i = 0; i < v.Length; i++)
        {
            v[i] += f[i] * dT / mass;
        }
        for (int i = 0; i < v.Length; i++)
        {
            vertices[i] += v[i] * dT;
        }
        vertices[N - 1] = new Vector3(0f, 0f, Size);
        vertices[(N - 1) * N + N - 1] = new Vector3(Size, 0f, Size);
        v[N - 1] = new Vector3(0f, 0f, 0f);
        v[(N - 1) * N + N - 1] = new Vector3(0f, 0f, 0f);
    }

    void GS()
    {
        float eps = 1e-3f;
        int n = X.Length;
        for (int i = 0; i < n; i++)
        {
            X[i] = 1f;
        }
        int I;
        for(I = 0; I < 30; I++)
        {
            float mx = 0f;
            for (int i = 0; i < n; i++)
            {
                float t = X[i];
                X[i] = B[i];
                float ii = 1;
                for (int j = 0; j < mat.data[i].Count; j++)
                {
                    if (mat.data[i][j].col != i)
                    {
                        X[i] -= mat.data[i][j].v * X[mat.data[i][j].col];
                    }
                    else
                    {
                        ii = mat.data[i][j].v;
                    }
                }
                X[i] /= ii;
                mx = Mathf.Max(mx, Mathf.Abs(X[i] - t));
            }
            if (mx < eps)
            {
                break;
            }
        }
        Debug.Log(I);
    }

    void AddSpringForce(int xj, int xi, float a, float b, ref Vector3 d, float c)
    {
        mat.Add(xi * 3 + 0, xi * 3 + 0, c * a * (-1f * (1f - b / d.magnitude) - b * Mathf.Pow(d.x, 2) / Mathf.Pow(d.magnitude, 3)));
        mat.Add(xi * 3 + 1, xi * 3 + 1, c * a * (-1f * (1f - b / d.magnitude) - b * Mathf.Pow(d.y, 2) / Mathf.Pow(d.magnitude, 3)));
        mat.Add(xi * 3 + 2, xi * 3 + 2, c * a * (-1f * (1f - b / d.magnitude) - b * Mathf.Pow(d.z, 2) / Mathf.Pow(d.magnitude, 3)));
        mat.Add(xi * 3 + 0, xj * 3 + 0, c * a * (1f * (1f - b / d.magnitude) + b * Mathf.Pow(d.x, 2) / Mathf.Pow(d.magnitude, 3)));
        mat.Add(xi * 3 + 1, xj * 3 + 1, c * a * (1f * (1f - b / d.magnitude) + b * Mathf.Pow(d.y, 2) / Mathf.Pow(d.magnitude, 3)));
        mat.Add(xi * 3 + 2, xj * 3 + 2, c * a * (1f * (1f - b / d.magnitude) + b * Mathf.Pow(d.z, 2) / Mathf.Pow(d.magnitude, 3)));
    }

    public void BackwardEuler(ref Vector3[] vertices)
    {
        tnow += dT;
        float lRestStructural = Size / (N - 1);
        float lRestShear = lRestStructural / Mathf.Sqrt(2);
        float lRestBend = lRestStructural * 2;

        float kStructural = kStructuralC / (1f / (N - 1));
        float kShear = kStructuralC / (1f / (N - 1) * Mathf.Sqrt(2f) / 2);
        float kBend = kBendC / (1f / (N - 1) * 2);
        mat = new MySparseMatrix(v.Length * 3);

        float mass = totMass / v.Length;
        float cc = -dT * dT / mass;
        for (int i = 0; i < N; i++)
        {
            for (int j = 0; j < N; j++)
            {
                f[i * N + j] = new Vector3(0f, 0f, 0f);
                if (i - 2 >= 0)
                {
                    Vector3 d = vertices[(i - 2) * N + j] - vertices[i * N + j];
                    f[i * N + j] += kBend * d.normalized * (d.magnitude - lRestBend);
                    AddSpringForce((i - 2) * N + j, i * N + j, kBend, lRestBend, ref d, cc);
                }
                if (j - 2 >= 0)
                {
                    Vector3 d = vertices[i * N + j - 2] - vertices[i * N + j];
                    f[i * N + j] += kBend * d.normalized * (d.magnitude - lRestBend);
                    AddSpringForce(i * N + j - 2, i * N + j, kBend, lRestBend, ref d, cc);
                }
                if (i + 2 < N)
                {
                    Vector3 d = vertices[(i + 2) * N + j] - vertices[i * N + j];
                    f[i * N + j] += kBend * d.normalized * (d.magnitude - lRestBend);
                    AddSpringForce((i + 2) * N + j, i * N + j, kBend, lRestBend, ref d, cc);
                }
                if (j + 2 < N)
                {
                    Vector3 d = vertices[i * N + j + 2] - vertices[i * N + j];
                    f[i * N + j] += kBend * d.normalized * (d.magnitude - lRestBend);
                    AddSpringForce(i * N + j + 2, i * N + j, kBend, lRestBend, ref d, cc);
                }

                if (i - 1 >= 0)
                {
                    Vector3 d = vertices[(i - 1) * N + j] - vertices[i * N + j];
                    f[i * N + j] += kStructural * d.normalized * (d.magnitude - lRestStructural);
                    AddSpringForce((i - 1) * N + j, i * N + j, kStructural, lRestStructural, ref d, cc);
                }
                if (j - 1 >= 0)
                {
                    Vector3 d = vertices[i * N + j - 1] - vertices[i * N + j];
                    f[i * N + j] += kStructural * d.normalized * (d.magnitude - lRestStructural);
                    AddSpringForce(i * N + j - 1, i * N + j, kStructural, lRestStructural, ref d, cc);
                }
                if (i + 1 < N)
                {
                    Vector3 d = vertices[(i + 1) * N + j] - vertices[i * N + j];
                    f[i * N + j] += kStructural * d.normalized * (d.magnitude - lRestStructural);
                    AddSpringForce((i + 1) * N + j, i * N + j, kStructural, lRestStructural, ref d, cc);
                }
                if (j + 1 < N)
                {
                    Vector3 d = vertices[i * N + j + 1] - vertices[i * N + j];
                    f[i * N + j] += kStructural * d.normalized * (d.magnitude - lRestStructural);
                    AddSpringForce(i * N + j + 1, i * N + j, kStructural, lRestStructural, ref d, cc);
                }
                if (i - 1 >= 0 && j - 1 >= 0)
                {
                    Vector3 d = vertices[N * N + (i - 1) * (N - 1) + j - 1] - vertices[i * N + j];
                    f[i * N + j] += kShear * d.normalized * (d.magnitude - lRestShear);
                    AddSpringForce(N * N + (i - 1) * (N - 1) + j - 1, i * N + j, kShear, lRestShear, ref d, cc);
                }
                if (i - 1 >= 0 && j < N - 1)
                {
                    Vector3 d = vertices[N * N + (i - 1) * (N - 1) + j] - vertices[i * N + j];
                    f[i * N + j] += kShear * d.normalized * (d.magnitude - lRestShear);
                    AddSpringForce(N * N + (i - 1) * (N - 1) + j, i * N + j, kShear, lRestShear, ref d, cc);
                }
                if (i < N - 1 && j - 1 >= 0)
                {
                    Vector3 d = vertices[N * N + i * (N - 1) + j - 1] - vertices[i * N + j];
                    f[i * N + j] += kShear * d.normalized * (d.magnitude - lRestShear);
                    AddSpringForce(N * N + i * (N - 1) + j - 1, i * N + j, kShear, lRestShear, ref d, cc);
                }
                if (i < N - 1 && j < N - 1)
                {
                    Vector3 d = vertices[N * N + i * (N - 1) + j] - vertices[i * N + j];
                    f[i * N + j] += kShear * d.normalized * (d.magnitude - lRestShear);
                    AddSpringForce(N * N + i * (N - 1) + j, i * N + j, kShear, lRestShear, ref d, cc);
                }
            }
        }
        for (int i = 0; i < N - 1; i++)
        {
            for (int j = 0; j < N - 1; j++)
            {
                f[N * N + i * (N - 1) + j] = new Vector3(0f, 0f, 0f);
                if (i - 2 >= 0)
                {
                    Vector3 d = vertices[N * N + (i - 2) * (N - 1) + j] - vertices[N * N + i * (N - 1) + j];
                    f[N * N + i * (N - 1) + j] += kBend * d.normalized * (d.magnitude - lRestBend);
                    AddSpringForce(N * N + (i - 2) * (N - 1) + j, N * N + i * (N - 1) + j, kBend, lRestBend, ref d, cc);
                }
                if (j - 2 >= 0)
                {
                    Vector3 d = vertices[N * N + i * (N - 1) + j - 2] - vertices[N * N + i * (N - 1) + j];
                    f[N * N + i * (N - 1) + j] += kBend * d.normalized * (d.magnitude - lRestBend);
                    AddSpringForce(N * N + i * (N - 1) + j - 2, N * N + i * (N - 1) + j, kBend, lRestBend, ref d, cc);
                }
                if (i + 2 < N - 1)
                {
                    Vector3 d = vertices[N * N + (i + 2) * (N - 1) + j] - vertices[N * N + i * (N - 1) + j];
                    f[N * N + i * (N - 1) + j] += kBend * d.normalized * (d.magnitude - lRestBend);
                    AddSpringForce(N * N + (i + 2) * (N - 1) + j, N * N + i * (N - 1) + j, kBend, lRestBend, ref d, cc);
                }
                if (j + 2 < N - 1)
                {
                    Vector3 d = vertices[N * N + i * (N - 1) + j + 2] - vertices[N * N + i * (N - 1) + j];
                    f[N * N + i * (N - 1) + j] += kBend * d.normalized * (d.magnitude - lRestBend);
                    AddSpringForce(N * N + i * (N - 1) + j + 2, N * N + i * (N - 1) + j, kBend, lRestBend, ref d, cc);
                }

                if (i - 1 >= 0)
                {
                    Vector3 d = vertices[N * N + (i - 1) * (N - 1) + j] - vertices[N * N + i * (N - 1) + j];
                    f[N * N + i * (N - 1) + j] += kStructural * d.normalized * (d.magnitude - lRestStructural);
                    AddSpringForce(N * N + (i - 1) * (N - 1) + j, N * N + i * (N - 1) + j, kStructural, lRestStructural, ref d, cc);
                }
                if (j - 1 >= 0)
                {
                    Vector3 d = vertices[N * N + i * (N - 1) + j - 1] - vertices[N * N + i * (N - 1) + j];
                    f[N * N + i * (N - 1) + j] += kStructural * d.normalized * (d.magnitude - lRestStructural);
                    AddSpringForce(N * N + i * (N - 1) + j - 1, N * N + i * (N - 1) + j, kStructural, lRestStructural, ref d, cc);
                }
                if (i + 1 < N - 1)
                {
                    Vector3 d = vertices[N * N + (i + 1) * (N - 1) + j] - vertices[N * N + i * (N - 1) + j];
                    f[N * N + i * (N - 1) + j] += kStructural * d.normalized * (d.magnitude - lRestStructural);
                    AddSpringForce(N * N + (i + 1) * (N - 1) + j, N * N + i * (N - 1) + j, kStructural, lRestStructural, ref d, cc);
                }
                if (j + 1 < N - 1)
                {
                    Vector3 d = vertices[N * N + i * (N - 1) + j + 1] - vertices[N * N + i * (N - 1) + j];
                    f[N * N + i * (N - 1) + j] += kStructural * d.normalized * (d.magnitude - lRestStructural);
                    AddSpringForce(N * N + i * (N - 1) + j + 1, N * N + i * (N - 1) + j, kStructural, lRestStructural, ref d, cc);
                }

                {
                    Vector3 d = vertices[i * N + j] - vertices[N * N + i * (N - 1) + j];
                    f[N * N + i * (N - 1) + j] += kShear * d.normalized * (d.magnitude - lRestShear);
                    AddSpringForce(i * N + j, N * N + i * (N - 1) + j, kShear, lRestShear, ref d, cc);
                }
                {
                    Vector3 d = vertices[i * N + j + 1] - vertices[N * N + i * (N - 1) + j];
                    f[N * N + i * (N - 1) + j] += kShear * d.normalized * (d.magnitude - lRestShear);
                    AddSpringForce(i * N + j + 1, N * N + i * (N - 1) + j, kShear, lRestShear, ref d, cc);
                }
                {
                    Vector3 d = vertices[(i + 1) * N + j] - vertices[N * N + i * (N - 1) + j];
                    f[N * N + i * (N - 1) + j] += kShear * d.normalized * (d.magnitude - lRestShear);
                    AddSpringForce((i + 1) * N + j, N * N + i * (N - 1) + j, kShear, lRestShear, ref d, cc);
                }
                {
                    Vector3 d = vertices[(i + 1) * N + j + 1] - vertices[N * N + i * (N - 1) + j];
                    f[N * N + i * (N - 1) + j] += kShear * d.normalized * (d.magnitude - lRestShear);
                    AddSpringForce((i + 1) * N + j + 1, N * N + i * (N - 1) + j, kShear, lRestShear, ref d, cc);
                }
            }
        }

        for (int i = 0; i < f.Length; i++)
        {
            f[i] += new Vector3(0f, -1f, 0f) * 9.8f * mass;
        }
        for (int i = 0; i < f.Length; i++)
        {
            f[i] += -new Vector3(v0[i * 3 + 0], v0[i * 3 + 1], v0[i * 3 + 2]) * kDamp;
        }

        for (int i = 0; i < f.Length; i++)
        {
            f[i] += new Vector3(0f, 0f, -1f) * Mathf.Pow(Mathf.Sin(tnow * tWind), 2f) * fWind;
        }

        mat.Multiply(ref v0, ref B);
        for (int i = 0; i < f.Length; i++)
        {
            B[i * 3 + 0] = dT / mass * (B[i * 3 + 0] * dT + f[i].x);
            B[i * 3 + 1] = dT / mass * (B[i * 3 + 1] * dT + f[i].y);
            B[i * 3 + 2] = dT / mass * (B[i * 3 + 2] * dT + f[i].z - dT * 2f * fWind * Mathf.Sin(tnow * tWind) * Mathf.Cos(tnow * tWind) * tWind);
        }
        for (int i = 0; i < v0.Length; i++)
        {
            mat.Add(i, i, 1 + kDamp / mass * dT);
        }

        B[(N - 1) * 3 + 0] = 0;
        B[(N - 1) * 3 + 1] = 0;
        B[(N - 1) * 3 + 2] = 0;
        B[((N - 1) * N + N - 1) * 3 + 0] = 0;
        B[((N - 1) * N + N - 1) * 3 + 1] = 0;
        B[((N - 1) * N + N - 1) * 3 + 2] = 0;


        mat.data[(N - 1) * 3 + 0] = new List<(int col, float v)>();
        mat.data[(N - 1) * 3 + 1] = new List<(int col, float v)>();
        mat.data[(N - 1) * 3 + 2] = new List<(int col, float v)>();
        mat.data[((N - 1) * N + N - 1) * 3 + 0] = new List<(int col, float v)>();
        mat.data[((N - 1) * N + N - 1) * 3 + 1] = new List<(int col, float v)>();
        mat.data[((N - 1) * N + N - 1) * 3 + 2] = new List<(int col, float v)>();

        mat.Add((N - 1) * 3 + 0, (N - 1) * 3 + 0, 1);
        mat.Add((N - 1) * 3 + 1, (N - 1) * 3 + 1, 1);
        mat.Add((N - 1) * 3 + 2, (N - 1) * 3 + 2, 1);
        mat.Add(((N - 1) * N + N - 1) * 3 + 0, (N - 1) * 3 + 0, 1);
        mat.Add(((N - 1) * N + N - 1) * 3 + 1, (N - 1) * 3 + 1, 1);
        mat.Add(((N - 1) * N + N - 1) * 3 + 2, (N - 1) * 3 + 2, 1);
        GS();
        for (int i = 0; i < v0.Length; i++)
        {
            v0[i] += X[i]; 
        }
        for (int i = 0; i < vertices.Length; i++)
        {
            vertices[i] += dT * new Vector3(v0[i * 3 + 0], v0[i * 3 + 1], v0[i * 3 + 2]);
        }
    }

    void Update()
    {
        Vector3[] vertices = gameObject.GetComponent<MeshFilter>().mesh.vertices;
        //ForwardEuler(ref vertices);
        BackwardEuler(ref vertices);
        List<Vector3> normals = new List<Vector3>();
        CalcNormal(ref vertices, ref normals);
        gameObject.GetComponent<MeshFilter>().mesh.vertices = vertices;
        gameObject.GetComponent<MeshFilter>().mesh.SetNormals(normals);
        gameObject.GetComponent<MeshFilter>().mesh.RecalculateBounds();
    }
}
