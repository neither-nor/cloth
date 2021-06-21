using System.Collections;
using System.Collections.Generic;
using UnityEngine;

public class CameraController : MonoBehaviour
{
    // Start is called before the first frame update
    void Start()
    {

    }


    public float sensitivity = 1f;
    public float speed = 3f;
    bool rotation = false;
    // Update is called once 
    void LateUpdate()
    {
        if (Input.GetMouseButtonDown(1))
        {
            rotation = true;
        }
        if (Input.GetMouseButtonUp(1))
        {
            rotation = false;
        }
        if (rotation)
        {
            float rotateHorizontal = Input.GetAxis("Mouse X");
            float rotateVertical = Input.GetAxis("Mouse Y");
            transform.RotateAround(transform.position, Vector3.up, rotateHorizontal * sensitivity); //use transform.Rotate(-transform.up * rotateHorizontal * sensitivity) instead if you dont want the camera to rotate around the player
            transform.RotateAround(transform.position, transform.right, -rotateVertical * sensitivity); // again, use transform.Rotate(transform.right * rotateVertical * sensitivity) if you don't want the camera to rotate around the player
        }
        float dz = Input.GetAxis("Horizontal");
        float dx = Input.GetAxis("Vertical");
        if (Mathf.Abs(dx) > 0 || Mathf.Abs(dz) > 0)
        {
            Vector3 translation = (dx * transform.forward + dz * transform.right).normalized;
            transform.position += Time.deltaTime * speed * translation;
        }
    }
}
