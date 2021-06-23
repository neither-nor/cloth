using System.Collections;
using System.Collections.Generic;
using UnityEngine;
using UnityEditor;

public class ClothManager : MonoBehaviour
{
    // Start is called before the first frame update
    public Mesh clothMesh;
    const float Size = 10f;
    [Range(3, 100)] public int N = 10;
    public float totMass = 100;
    Vector3[] f;
    Vector3[] v;
    float[] v0;
    float[] dv0;
    float[] f0;
    float[] B;
    float[] X;
    float dT = 0.005f;
    List<int> triangles;
    List<Vector2> uvs;
    List<Vector3> verticesRest;
    public float kStructuralC = 10f;
    public float kBendC = 1f;
    public float kDamp = 5f;
    public float fWind = 1f;
    public float tWind = 1f;
    public bool fixStrain;

    float lRestStructural;
    float lRestShear;
    float lRestBend;

    float kStructural;
    float kShear;
    float kBend;

    float tnow = 0;
    MySparseMatrix mat;
    void Start()
    {
        f = new Vector3[N * N];
        v = new Vector3[N * N];
        f0 = new float[f.Length * 3];
        v0 = new float[v.Length * 3];
        dv0 = new float[v.Length * 3];
        B = new float[v.Length * 3];
        X = new float[v.Length * 3];


        lRestStructural = Size / (N - 1);
        lRestShear = lRestStructural * Mathf.Sqrt(2);
        lRestBend = lRestStructural * 2;

        kStructural = kStructuralC / (1f / (N - 1));
        kShear = kStructuralC / (1f / (N - 1) * Mathf.Sqrt(2f));
        kBend = kBendC / (1f / (N - 1) * 2);

        BuildClothMesh();
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
    }

    public void BuildClothMesh()
    {
        MeshFilter meshFilter = gameObject.GetComponent<MeshFilter>();

        clothMesh = new Mesh();
        verticesRest = new List<Vector3>();
        uvs = new List<Vector2>();
        List<Vector3> normals = new List<Vector3>();
        triangles = new List<int>();
        for (int i = 0; i < N; i++)
        {
            for (int j = 0; j < N; j++)
            {
                Vector2 uv = new Vector2(1f * i / (N - 1), 1f * j / (N - 1));
                uvs.Add(uv);
                Vector3 v = new Vector3(uv.x * Size,
                                        0f,
                                        uv.y * Size);
                verticesRest.Add(v);
            }
        }
        for (int i = 0; i < N - 1; i++)
        {
            for (int j = 0; j < N - 1; j++)
            {
                triangles.Add(i * N + j);
                triangles.Add(i * N + j + 1);
                triangles.Add(i * N + j + 1 + N);

                triangles.Add(i * N + j);
                triangles.Add(i * N + j + 1 + N);
                triangles.Add(i * N + j + N);

            }
        }

        //CalcNormal(ref vertices, ref normals);

        clothMesh.SetVertices(verticesRest);
        clothMesh.SetUVs(0, uvs);
        clothMesh.SetTriangles(triangles, 0);
        clothMesh.SetNormals(normals);
        meshFilter.mesh = clothMesh;
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
                    Vector3 d = vertices[(i - 1) * N + j - 1] - vertices[i * N + j];
                    f[i * N + j] += kShear * d.normalized * (d.magnitude - lRestShear);
                }
                if (i - 1 >= 0 && j + 1 < N)
                {
                    Vector3 d = vertices[(i - 1) * N + j + 1] - vertices[i * N + j];
                    f[i * N + j] += kShear * d.normalized * (d.magnitude - lRestShear);
                }
                if (i + 1 < N && j - 1 >= 0)
                {
                    Vector3 d = vertices[(i + 1) * N + j - 1] - vertices[i * N + j];
                    f[i * N + j] += kShear * d.normalized * (d.magnitude - lRestShear);
                }
                if (i + 1 < N && j + 1 < N)
                {
                    Vector3 d = vertices[(i + 1) * N + j + 1] - vertices[i * N + j];
                    f[i * N + j] += kShear * d.normalized * (d.magnitude - lRestShear);
                }
            }
        }
        float mass = totMass / v.Length;
        for (int i = 0; i < f.Length; i++)
        {
            f[i] += new Vector3(0f, -1f, 0f) * 9.8f * mass;
        }
        for (int i = 0; i < f.Length; i++)
        {
            f[i] += -v[i] * kDamp * Size * Size / (N * N);
        }

        tnow += dT;
        for (int i = 0; i < f.Length; i++)
        {
            f[i] += new Vector3(0f, 0f, -1f) * Mathf.Pow(Mathf.Sin(tnow * tWind), 2f) * fWind * Size * Size / (N * N);
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
        for (I = 0; I < 30; I++)
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

    Vector3 V0(int x) {
        return new Vector3(v0[x * 3 + 0], v0[x * 3 + 1], v0[x * 3 + 2]);
    }

    void fvadd(ref float[] fv, int x, Vector3 v) {
        fv[x * 3 + 0] += v.x;
        fv[x * 3 + 1] += v.y;
        fv[x * 3 + 2] += v.z;
    }

    void fvset(ref float[] fv, int x, Vector3 v)
    {
        fv[x * 3 + 0] = v.x;
        fv[x * 3 + 1] = v.y;
        fv[x * 3 + 2] = v.z;
    }

    bool FixEdgeStrain(int x, int y, ref Vector3[] vertices) {
        float lRest = (verticesRest[x] - verticesRest[y]).magnitude;
        Vector3 d = (vertices[x] + dT * V0(x)) - (vertices[y] + dT * V0(y));
        float lNow = d.magnitude;
        if (lNow > lRest * 1.1f)
        {
            Vector3 cv = (lNow - lRest * 1.1f) / dT / 2f * d.normalized / 2f;
            fvadd(ref dv0, x, -cv);
            fvadd(ref dv0, y, cv);
            return true;
        }
        if (lNow < lRest * .9f)
        {
            Vector3 cv = (lNow - lRest * .9f) / dT / 2f * d.normalized / 2f;
            fvadd(ref dv0, x, -cv);
            fvadd(ref dv0, y, cv);
            return true;
        }
        return false;
    }

    void FixStrain(ref Vector3[] vertices)
    {
        bool flag = true;
        for (int i = 0; i < dv0.Length; i++)
        {
            dv0[i] = 0;
        }
        int I;
        for(I = 0; I < 5; I++)
        {
            flag = false;
            for (int i = 0; i < triangles.Count; i += 3)
            {
                flag = flag || FixEdgeStrain(triangles[i], triangles[i + 1], ref vertices);
            }
            if (flag)
            {
                for (int i = 0; i < v0.Length; i++)
                {
                    v0[i] += dv0[i];
                    dv0[i] = 0;
                }
            }
            else
            {
                break;
            }
        }
        Debug.LogFormat("FixStrain{0}", I);
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

    public void SpringMassBackwardEuler(ref Vector3[] vertices)
    {
        tnow += dT;


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
                    Vector3 d = vertices[(i - 1) * N + j - 1] - vertices[i * N + j];
                    f[i * N + j] += kShear * d.normalized * (d.magnitude - lRestShear);
                    AddSpringForce((i - 1) * N + j - 1, i * N + j, kShear, lRestShear, ref d, cc);
                }
                if (i - 1 >= 0 && j + 1 < N)
                {
                    Vector3 d = vertices[(i - 1) * N + j + 1] - vertices[i * N + j];
                    f[i * N + j] += kShear * d.normalized * (d.magnitude - lRestShear);
                    AddSpringForce((i - 1) * N + j + 1, i * N + j, kShear, lRestShear, ref d, cc);
                }
                if (i + 1 < N && j - 1 >= 0)
                {
                    Vector3 d = vertices[(i + 1) * N + j - 1] - vertices[i * N + j];
                    f[i * N + j] += kShear * d.normalized * (d.magnitude - lRestShear);
                    AddSpringForce((i + 1) * N + j - 1, i * N + j, kShear, lRestShear, ref d, cc);
                }
                if (i + 1 < N && j + 1 < N)
                {
                    Vector3 d = vertices[(i + 1) * N + j + 1] - vertices[i * N + j];
                    f[i * N + j] += kShear * d.normalized * (d.magnitude - lRestShear);
                    AddSpringForce((i + 1) * N + j + 1, i * N + j, kShear, lRestShear, ref d, cc);
                }
            }
        }

        for (int i = 0; i < f.Length; i++)
        {
            f[i] += new Vector3(0f, -1f, 0f) * 9.8f * mass;
        }
        for (int i = 0; i < f.Length; i++)
        {
            f[i] += -new Vector3(v0[i * 3 + 0], v0[i * 3 + 1], v0[i * 3 + 2]) * kDamp * Size * Size / (N * N);
        }

        for (int i = 0; i < f.Length; i++)
        {
            f[i] += new Vector3(0f, 0f, -1f) * Mathf.Pow(Mathf.Sin(tnow * tWind), 2f) * fWind * Size * Size / (N * N);
        }

        mat.Multiply(ref v0, ref B);
        for (int i = 0; i < f.Length; i++)
        {
            B[i * 3 + 0] = dT / mass * (B[i * 3 + 0] * dT + f[i].x);
            B[i * 3 + 1] = dT / mass * (B[i * 3 + 1] * dT + f[i].y);
            B[i * 3 + 2] = dT / mass * (B[i * 3 + 2] * dT + f[i].z - dT * 2f * fWind * Mathf.Sin(tnow * tWind) * Mathf.Cos(tnow * tWind) * tWind * Size * Size / (N * N));
        }
        for (int i = 0; i < v0.Length; i++)
        {
            mat.Add(i, i, 1 + kDamp / mass * dT * Size * Size / (N * N));
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
        if (fixStrain) {
            FixStrain(ref vertices);
            vertices[N - 1] = new Vector3(0f, 0f, Size);
            vertices[(N - 1) * N + N - 1] = new Vector3(Size, 0f, Size);
            fvset(ref v0, N - 1, new Vector3(0f, 0f, 0f));
            fvset(ref v0, (N - 1) * N + N - 1, new Vector3(0f, 0f, 0f));
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
        SpringMassBackwardEuler(ref vertices);
        List<Vector3> normals = new List<Vector3>();
        CalcNormal(ref vertices, ref normals);
        gameObject.GetComponent<MeshFilter>().mesh.vertices = vertices;
        gameObject.GetComponent<MeshFilter>().mesh.SetNormals(normals);
        gameObject.GetComponent<MeshFilter>().mesh.RecalculateBounds();
    }
}
