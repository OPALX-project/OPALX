import os, sys, argparse, math, base64, zlib, struct

if sys.version_info < (3,0):
    range = xrange

vertices = []
numVertices = [74, 74, 74, 74, 74, 74, 74, 74, 74, 74, 74, 74, 74, 80, 74, 74, 74, 74, 80, 74, 74, 74, 74, 74, 104, 74, 74, 74]
vertices_base64 = "eJztnflfjdv7/yuVJBWRKGSmRJKhxL5pPMYM5xwZk8zzPM/zPM9T5hAiFMq9zZkrlcxFFCJqIwrftbled7rOt8fnD3jXb+fpPqu1rvW61rPV3pWOTsGPpotb+b/0laTf/5Ws4ry8zexlS+c+Uy0eW3JGRdfTKnAbh86tP5d+qvpUcsm44BdXFH64Xbf0KacfqX7++rircBffu7tHzk1S2Vi6dQ5fHa/wSvTf+Hdw3V+ziFEdovHAv9Dns6bPD76M5of5ghfvWbZRrf2z3bA+8EW/n5f587Qe2YaNn/d7PfJhNp+Kv+cvu7H5N/n933JVtt4Dv/9/uTitD9zq9+eTi5sUrGe53/OTT7P603rk2UsDrl52vuGmzIee70jrA69K4yfS+sBP0nySaX3gXjT/UrQ+cHv6bwf6d/DSv9cjB9F44FV+r0cuS58fPI3mZ0HzBd+Vfn3d2gMPIrE+JSe0fhv2/HnKQy02/g3K31k2n2K0Hx3Y/Ktj/9h6N9F+29L6wE0pHy1NCtbTnPJUd1zB+uv8/z5m5fed2rNhv4CKLyPBq9A4j8f8Xje4PX3enbRu8Is0zyO0bnA/WteXcr/XDe5BdbCifwevT7lcReOB96JclqDPDz7k9zpVpWi+4A03fs/95v9O6Tvws7/XI9uy51dSLuuz8XdSLq+x+Xz4vR65D5u/Be1rW8uC611GOWhC6wM3pNwMNylYTxPK2bJxBeuPXPpqAm9cv/fBDbw6PT+Z1gfuTOOPovWB36X5zKH1gQ+m+b+g9YH/S7k0on9X9pFyOZfGU9ZLudShzw9+m3JZguYLblTppqXlF43Sd+DNKH+12POdKH+ubPz+lL/7bD73KX/j2Py/d/y9f4GWBdc7lfa7M60PPLfe73wcNilYT0PKk+n4gvXXYR/+L2OqPdr8X9+Bc9+Bc9+Bc9+Bc9+Bc9+Bc9+Bc9+Bc9+Bc9+Bc9+Bc9+Bc9+Bc9+Bc9+Bc9+Bc9+Bc9+x/VJ8p8yH+Q6c+w6c+w6c+w6c+w6c+w6c+w6c+w6c+07JCfMdOPcdOPcdOPcdOPcdOPcdOPcdOPddwf3K/6jbc+XFT9H/9R049x049x049x049x049x049x049x049x049x049x049x049x049x049x049x049x049x3bL8V34Nx34Nx34Nx34Nx34Nx3yj5SLuE7Zb3Md+Dcd+Dcd+Dcd+Dcd+Dcd+Dcd+Dcd+Dcd+Dcd8o+Mt8V3K/8D+2XnzHP/+s7cO47cO47cO47cO47cO47cO47cO47cO47cO47cO47cO47cO47cO47cO47cO47cO47cO47tl+K75T5MN+Bc9+Bc9+Bc9+Bc9+Bc9+Bc9+Bc9+Bc98pOWG+A+e+A+e+A+e+A+e+A+e+A+e+A+e+K7hf+R+dIyrPf1u61X98B859B859B859B859B859B859B859B859B859B859B859B859B859B859B859B859B859x/ZL8R049x049x049x049x04952yj8x3ynqZ78C578C578C578C578C578C578C578C578C575R9ZL4ruF/5H5203CW/7zjnvgPnvgPnvgPnvgPnvgPnvgPnvgPnvgPnvgPnvgPnvgPnvgPnvgPnvgPnvgPnvgPnvmP7pfhOmQ/zHTj3HTj3HTj3HTj3HTj3HTj3HTj3HTj3nZIT5jtw7jtw7jtw7jtw7jtw7jtw7jtw7ruC+5X/4dRIfDT/r+/Aue/Aue/Aue/Aue/Aue/Aue/Aue/Aue/Aue/Aue/Aue/Aue/Aue/Aue/Aue/Aue/Aue/Yfim+A+e+A+e+A+e+A+e+A+e+U/aR+U5ZL/MdOPcdOPcdOPcdOPcdOPcdOPcdOPcdOPcdOPedso/MdwX3q8h3Rb4r8l2R74p8V+S7/y3fPb+26dwGq9b/8R049x049x049x049x049x049x049x049x049x049x049x049x049x049x049x049x049x3bL8V3ynyY78C578C578C578C578C578C578C578C575ScMN+Bc9+Bc9+Bc9+Bc9+Bc9+Bc9+Bc98V3K/8j0utjwx+WiG/7+A7cO47cO47cO47cO47cO47cO47cO47cO47cO47cO47cO47cO47cO47cO47cO47cO47cO47tl+K78C578C578C578C578C575R9ZL5T1st8B859B859B859B859B859B859B859B859B859p+wj813B/SryXZHvinxX5Lsi3xX57n/Ld452By1iavzXd+Dcd+Dcd+Dcd+Dcd+Dcd+Dcd+Dcd+Dcd+Dcd+Dcd+Dcd+Dcd+Dcd+Dcd+Dcd+Dcd+Dcd2y/FN8p82G+A+e+A+e+A+e+A+e+A+e+A+e+A+e+A+e+U3LCfAfOfQfO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"
triangles = [[ 0, 2, 1, 1, 2, 38, 2, 39, 38, 0, 3, 2, 2, 3, 39, 3, 40, 39, 0, 4, 3, 3, 4, 40, 4, 41, 40, 0, 5, 4, 4, 5, 41, 5, 42, 41, 0, 6, 5, 5, 6, 42, 6, 43, 42, 0, 7, 6, 6, 7, 43, 7, 44, 43, 0, 8, 7, 7, 8, 44, 8, 45, 44, 0, 9, 8, 8, 9, 45, 9, 46, 45, 0, 10, 9, 9, 10, 46, 10, 47, 46, 0, 11, 10, 10, 11, 47, 11, 48, 47, 0, 12, 11, 11, 12, 48, 12, 49, 48, 0, 13, 12, 12, 13, 49, 13, 50, 49, 0, 14, 13, 13, 14, 50, 14, 51, 50, 0, 15, 14, 14, 15, 51, 15, 52, 51, 0, 16, 15, 15, 16, 52, 16, 53, 52, 0, 17, 16, 16, 17, 53, 17, 54, 53, 0, 18, 17, 17, 18, 54, 18, 55, 54, 0, 19, 18, 18, 19, 55, 19, 56, 55, 0, 20, 19, 19, 20, 56, 20, 57, 56, 0, 21, 20, 20, 21, 57, 21, 58, 57, 0, 22, 21, 21, 22, 58, 22, 59, 58, 0, 23, 22, 22, 23, 59, 23, 60, 59, 0, 24, 23, 23, 24, 60, 24, 61, 60, 0, 25, 24, 24, 25, 61, 25, 62, 61, 0, 26, 25, 25, 26, 62, 26, 63, 62, 0, 27, 26, 26, 27, 63, 27, 64, 63, 0, 28, 27, 27, 28, 64, 28, 65, 64, 0, 29, 28, 28, 29, 65, 29, 66, 65, 0, 30, 29, 29, 30, 66, 30, 67, 66, 0, 31, 30, 30, 31, 67, 31, 68, 67, 0, 32, 31, 31, 32, 68, 32, 69, 68, 0, 33, 32, 32, 33, 69, 33, 70, 69, 0, 34, 33, 33, 34, 70, 34, 71, 70, 0, 35, 34, 34, 35, 71, 35, 72, 71, 0, 36, 35, 35, 36, 72, 36, 73, 72, 0, 1, 36, 36, 1, 73, 1, 38, 73, 37, 38, 39, 37, 39, 40, 37, 40, 41, 37, 41, 42, 37, 42, 43, 37, 43, 44, 37, 44, 45, 37, 45, 46, 37, 46, 47, 37, 47, 48, 37, 48, 49, 37, 49, 50, 37, 50, 51, 37, 51, 52, 37, 52, 53, 37, 53, 54, 37, 54, 55, 37, 55, 56, 37, 56, 57, 37, 57, 58, 37, 58, 59, 37, 59, 60, 37, 60, 61, 37, 61, 62, 37, 62, 63, 37, 63, 64, 37, 64, 65, 37, 65, 66, 37, 66, 67, 37, 67, 68, 37, 68, 69, 37, 69, 70, 37, 70, 71, 37, 71, 72, 37, 72, 73, 37, 73, 38], [ 0, 2, 1, 1, 2, 38, 2, 39, 38, 0, 3, 2, 2, 3, 39, 3, 40, 39, 0, 4, 3, 3, 4, 40, 4, 41, 40, 0, 5, 4, 4, 5, 41, 5, 42, 41, 0, 6, 5, 5, 6, 42, 6, 43, 42, 0, 7, 6, 6, 7, 43, 7, 44, 43, 0, 8, 7, 7, 8, 44, 8, 45, 44, 0, 9, 8, 8, 9, 45, 9, 46, 45, 0, 10, 9, 9, 10, 46, 10, 47, 46, 0, 11, 10, 10, 11, 47, 11, 48, 47, 0, 12, 11, 11, 12, 48, 12, 49, 48, 0, 13, 12, 12, 13, 49, 13, 50, 49, 0, 14, 13, 13, 14, 50, 14, 51, 50, 0, 15, 14, 14, 15, 51, 15, 52, 51, 0, 16, 15, 15, 16, 52, 16, 53, 52, 0, 17, 16, 16, 17, 53, 17, 54, 53, 0, 18, 17, 17, 18, 54, 18, 55, 54, 0, 19, 18, 18, 19, 55, 19, 56, 55, 0, 20, 19, 19, 20, 56, 20, 57, 56, 0, 21, 20, 20, 21, 57, 21, 58, 57, 0, 22, 21, 21, 22, 58, 22, 59, 58, 0, 23, 22, 22, 23, 59, 23, 60, 59, 0, 24, 23, 23, 24, 60, 24, 61, 60, 0, 25, 24, 24, 25, 61, 25, 62, 61, 0, 26, 25, 25, 26, 62, 26, 63, 62, 0, 27, 26, 26, 27, 63, 27, 64, 63, 0, 28, 27, 27, 28, 64, 28, 65, 64, 0, 29, 28, 28, 29, 65, 29, 66, 65, 0, 30, 29, 29, 30, 66, 30, 67, 66, 0, 31, 30, 30, 31, 67, 31, 68, 67, 0, 32, 31, 31, 32, 68, 32, 69, 68, 0, 33, 32, 32, 33, 69, 33, 70, 69, 0, 34, 33, 33, 34, 70, 34, 71, 70, 0, 35, 34, 34, 35, 71, 35, 72, 71, 0, 36, 35, 35, 36, 72, 36, 73, 72, 0, 1, 36, 36, 1, 73, 1, 38, 73, 37, 38, 39, 37, 39, 40, 37, 40, 41, 37, 41, 42, 37, 42, 43, 37, 43, 44, 37, 44, 45, 37, 45, 46, 37, 46, 47, 37, 47, 48, 37, 48, 49, 37, 49, 50, 37, 50, 51, 37, 51, 52, 37, 52, 53, 37, 53, 54, 37, 54, 55, 37, 55, 56, 37, 56, 57, 37, 57, 58, 37, 58, 59, 37, 59, 60, 37, 60, 61, 37, 61, 62, 37, 62, 63, 37, 63, 64, 37, 64, 65, 37, 65, 66, 37, 66, 67, 37, 67, 68, 37, 68, 69, 37, 69, 70, 37, 70, 71, 37, 71, 72, 37, 72, 73, 37, 73, 38], [ 0, 2, 1, 1, 2, 38, 2, 39, 38, 0, 3, 2, 2, 3, 39, 3, 40, 39, 0, 4, 3, 3, 4, 40, 4, 41, 40, 0, 5, 4, 4, 5, 41, 5, 42, 41, 0, 6, 5, 5, 6, 42, 6, 43, 42, 0, 7, 6, 6, 7, 43, 7, 44, 43, 0, 8, 7, 7, 8, 44, 8, 45, 44, 0, 9, 8, 8, 9, 45, 9, 46, 45, 0, 10, 9, 9, 10, 46, 10, 47, 46, 0, 11, 10, 10, 11, 47, 11, 48, 47, 0, 12, 11, 11, 12, 48, 12, 49, 48, 0, 13, 12, 12, 13, 49, 13, 50, 49, 0, 14, 13, 13, 14, 50, 14, 51, 50, 0, 15, 14, 14, 15, 51, 15, 52, 51, 0, 16, 15, 15, 16, 52, 16, 53, 52, 0, 17, 16, 16, 17, 53, 17, 54, 53, 0, 18, 17, 17, 18, 54, 18, 55, 54, 0, 19, 18, 18, 19, 55, 19, 56, 55, 0, 20, 19, 19, 20, 56, 20, 57, 56, 0, 21, 20, 20, 21, 57, 21, 58, 57, 0, 22, 21, 21, 22, 58, 22, 59, 58, 0, 23, 22, 22, 23, 59, 23, 60, 59, 0, 24, 23, 23, 24, 60, 24, 61, 60, 0, 25, 24, 24, 25, 61, 25, 62, 61, 0, 26, 25, 25, 26, 62, 26, 63, 62, 0, 27, 26, 26, 27, 63, 27, 64, 63, 0, 28, 27, 27, 28, 64, 28, 65, 64, 0, 29, 28, 28, 29, 65, 29, 66, 65, 0, 30, 29, 29, 30, 66, 30, 67, 66, 0, 31, 30, 30, 31, 67, 31, 68, 67, 0, 32, 31, 31, 32, 68, 32, 69, 68, 0, 33, 32, 32, 33, 69, 33, 70, 69, 0, 34, 33, 33, 34, 70, 34, 71, 70, 0, 35, 34, 34, 35, 71, 35, 72, 71, 0, 36, 35, 35, 36, 72, 36, 73, 72, 0, 1, 36, 36, 1, 73, 1, 38, 73, 37, 38, 39, 37, 39, 40, 37, 40, 41, 37, 41, 42, 37, 42, 43, 37, 43, 44, 37, 44, 45, 37, 45, 46, 37, 46, 47, 37, 47, 48, 37, 48, 49, 37, 49, 50, 37, 50, 51, 37, 51, 52, 37, 52, 53, 37, 53, 54, 37, 54, 55, 37, 55, 56, 37, 56, 57, 37, 57, 58, 37, 58, 59, 37, 59, 60, 37, 60, 61, 37, 61, 62, 37, 62, 63, 37, 63, 64, 37, 64, 65, 37, 65, 66, 37, 66, 67, 37, 67, 68, 37, 68, 69, 37, 69, 70, 37, 70, 71, 37, 71, 72, 37, 72, 73, 37, 73, 38], [ 0, 2, 1, 1, 2, 38, 2, 39, 38, 0, 3, 2, 2, 3, 39, 3, 40, 39, 0, 4, 3, 3, 4, 40, 4, 41, 40, 0, 5, 4, 4, 5, 41, 5, 42, 41, 0, 6, 5, 5, 6, 42, 6, 43, 42, 0, 7, 6, 6, 7, 43, 7, 44, 43, 0, 8, 7, 7, 8, 44, 8, 45, 44, 0, 9, 8, 8, 9, 45, 9, 46, 45, 0, 10, 9, 9, 10, 46, 10, 47, 46, 0, 11, 10, 10, 11, 47, 11, 48, 47, 0, 12, 11, 11, 12, 48, 12, 49, 48, 0, 13, 12, 12, 13, 49, 13, 50, 49, 0, 14, 13, 13, 14, 50, 14, 51, 50, 0, 15, 14, 14, 15, 51, 15, 52, 51, 0, 16, 15, 15, 16, 52, 16, 53, 52, 0, 17, 16, 16, 17, 53, 17, 54, 53, 0, 18, 17, 17, 18, 54, 18, 55, 54, 0, 19, 18, 18, 19, 55, 19, 56, 55, 0, 20, 19, 19, 20, 56, 20, 57, 56, 0, 21, 20, 20, 21, 57, 21, 58, 57, 0, 22, 21, 21, 22, 58, 22, 59, 58, 0, 23, 22, 22, 23, 59, 23, 60, 59, 0, 24, 23, 23, 24, 60, 24, 61, 60, 0, 25, 24, 24, 25, 61, 25, 62, 61, 0, 26, 25, 25, 26, 62, 26, 63, 62, 0, 27, 26, 26, 27, 63, 27, 64, 63, 0, 28, 27, 27, 28, 64, 28, 65, 64, 0, 29, 28, 28, 29, 65, 29, 66, 65, 0, 30, 29, 29, 30, 66, 30, 67, 66, 0, 31, 30, 30, 31, 67, 31, 68, 67, 0, 32, 31, 31, 32, 68, 32, 69, 68, 0, 33, 32, 32, 33, 69, 33, 70, 69, 0, 34, 33, 33, 34, 70, 34, 71, 70, 0, 35, 34, 34, 35, 71, 35, 72, 71, 0, 36, 35, 35, 36, 72, 36, 73, 72, 0, 1, 36, 36, 1, 73, 1, 38, 73, 37, 38, 39, 37, 39, 40, 37, 40, 41, 37, 41, 42, 37, 42, 43, 37, 43, 44, 37, 44, 45, 37, 45, 46, 37, 46, 47, 37, 47, 48, 37, 48, 49, 37, 49, 50, 37, 50, 51, 37, 51, 52, 37, 52, 53, 37, 53, 54, 37, 54, 55, 37, 55, 56, 37, 56, 57, 37, 57, 58, 37, 58, 59, 37, 59, 60, 37, 60, 61, 37, 61, 62, 37, 62, 63, 37, 63, 64, 37, 64, 65, 37, 65, 66, 37, 66, 67, 37, 67, 68, 37, 68, 69, 37, 69, 70, 37, 70, 71, 37, 71, 72, 37, 72, 73, 37, 73, 38], [ 0, 2, 1, 1, 2, 38, 2, 39, 38, 0, 3, 2, 2, 3, 39, 3, 40, 39, 0, 4, 3, 3, 4, 40, 4, 41, 40, 0, 5, 4, 4, 5, 41, 5, 42, 41, 0, 6, 5, 5, 6, 42, 6, 43, 42, 0, 7, 6, 6, 7, 43, 7, 44, 43, 0, 8, 7, 7, 8, 44, 8, 45, 44, 0, 9, 8, 8, 9, 45, 9, 46, 45, 0, 10, 9, 9, 10, 46, 10, 47, 46, 0, 11, 10, 10, 11, 47, 11, 48, 47, 0, 12, 11, 11, 12, 48, 12, 49, 48, 0, 13, 12, 12, 13, 49, 13, 50, 49, 0, 14, 13, 13, 14, 50, 14, 51, 50, 0, 15, 14, 14, 15, 51, 15, 52, 51, 0, 16, 15, 15, 16, 52, 16, 53, 52, 0, 17, 16, 16, 17, 53, 17, 54, 53, 0, 18, 17, 17, 18, 54, 18, 55, 54, 0, 19, 18, 18, 19, 55, 19, 56, 55, 0, 20, 19, 19, 20, 56, 20, 57, 56, 0, 21, 20, 20, 21, 57, 21, 58, 57, 0, 22, 21, 21, 22, 58, 22, 59, 58, 0, 23, 22, 22, 23, 59, 23, 60, 59, 0, 24, 23, 23, 24, 60, 24, 61, 60, 0, 25, 24, 24, 25, 61, 25, 62, 61, 0, 26, 25, 25, 26, 62, 26, 63, 62, 0, 27, 26, 26, 27, 63, 27, 64, 63, 0, 28, 27, 27, 28, 64, 28, 65, 64, 0, 29, 28, 28, 29, 65, 29, 66, 65, 0, 30, 29, 29, 30, 66, 30, 67, 66, 0, 31, 30, 30, 31, 67, 31, 68, 67, 0, 32, 31, 31, 32, 68, 32, 69, 68, 0, 33, 32, 32, 33, 69, 33, 70, 69, 0, 34, 33, 33, 34, 70, 34, 71, 70, 0, 35, 34, 34, 35, 71, 35, 72, 71, 0, 36, 35, 35, 36, 72, 36, 73, 72, 0, 1, 36, 36, 1, 73, 1, 38, 73, 37, 38, 39, 37, 39, 40, 37, 40, 41, 37, 41, 42, 37, 42, 43, 37, 43, 44, 37, 44, 45, 37, 45, 46, 37, 46, 47, 37, 47, 48, 37, 48, 49, 37, 49, 50, 37, 50, 51, 37, 51, 52, 37, 52, 53, 37, 53, 54, 37, 54, 55, 37, 55, 56, 37, 56, 57, 37, 57, 58, 37, 58, 59, 37, 59, 60, 37, 60, 61, 37, 61, 62, 37, 62, 63, 37, 63, 64, 37, 64, 65, 37, 65, 66, 37, 66, 67, 37, 67, 68, 37, 68, 69, 37, 69, 70, 37, 70, 71, 37, 71, 72, 37, 72, 73, 37, 73, 38], [ 0, 2, 1, 1, 2, 38, 2, 39, 38, 0, 3, 2, 2, 3, 39, 3, 40, 39, 0, 4, 3, 3, 4, 40, 4, 41, 40, 0, 5, 4, 4, 5, 41, 5, 42, 41, 0, 6, 5, 5, 6, 42, 6, 43, 42, 0, 7, 6, 6, 7, 43, 7, 44, 43, 0, 8, 7, 7, 8, 44, 8, 45, 44, 0, 9, 8, 8, 9, 45, 9, 46, 45, 0, 10, 9, 9, 10, 46, 10, 47, 46, 0, 11, 10, 10, 11, 47, 11, 48, 47, 0, 12, 11, 11, 12, 48, 12, 49, 48, 0, 13, 12, 12, 13, 49, 13, 50, 49, 0, 14, 13, 13, 14, 50, 14, 51, 50, 0, 15, 14, 14, 15, 51, 15, 52, 51, 0, 16, 15, 15, 16, 52, 16, 53, 52, 0, 17, 16, 16, 17, 53, 17, 54, 53, 0, 18, 17, 17, 18, 54, 18, 55, 54, 0, 19, 18, 18, 19, 55, 19, 56, 55, 0, 20, 19, 19, 20, 56, 20, 57, 56, 0, 21, 20, 20, 21, 57, 21, 58, 57, 0, 22, 21, 21, 22, 58, 22, 59, 58, 0, 23, 22, 22, 23, 59, 23, 60, 59, 0, 24, 23, 23, 24, 60, 24, 61, 60, 0, 25, 24, 24, 25, 61, 25, 62, 61, 0, 26, 25, 25, 26, 62, 26, 63, 62, 0, 27, 26, 26, 27, 63, 27, 64, 63, 0, 28, 27, 27, 28, 64, 28, 65, 64, 0, 29, 28, 28, 29, 65, 29, 66, 65, 0, 30, 29, 29, 30, 66, 30, 67, 66, 0, 31, 30, 30, 31, 67, 31, 68, 67, 0, 32, 31, 31, 32, 68, 32, 69, 68, 0, 33, 32, 32, 33, 69, 33, 70, 69, 0, 34, 33, 33, 34, 70, 34, 71, 70, 0, 35, 34, 34, 35, 71, 35, 72, 71, 0, 36, 35, 35, 36, 72, 36, 73, 72, 0, 1, 36, 36, 1, 73, 1, 38, 73, 37, 38, 39, 37, 39, 40, 37, 40, 41, 37, 41, 42, 37, 42, 43, 37, 43, 44, 37, 44, 45, 37, 45, 46, 37, 46, 47, 37, 47, 48, 37, 48, 49, 37, 49, 50, 37, 50, 51, 37, 51, 52, 37, 52, 53, 37, 53, 54, 37, 54, 55, 37, 55, 56, 37, 56, 57, 37, 57, 58, 37, 58, 59, 37, 59, 60, 37, 60, 61, 37, 61, 62, 37, 62, 63, 37, 63, 64, 37, 64, 65, 37, 65, 66, 37, 66, 67, 37, 67, 68, 37, 68, 69, 37, 69, 70, 37, 70, 71, 37, 71, 72, 37, 72, 73, 37, 73, 38], [ 0, 2, 1, 1, 2, 38, 2, 39, 38, 0, 3, 2, 2, 3, 39, 3, 40, 39, 0, 4, 3, 3, 4, 40, 4, 41, 40, 0, 5, 4, 4, 5, 41, 5, 42, 41, 0, 6, 5, 5, 6, 42, 6, 43, 42, 0, 7, 6, 6, 7, 43, 7, 44, 43, 0, 8, 7, 7, 8, 44, 8, 45, 44, 0, 9, 8, 8, 9, 45, 9, 46, 45, 0, 10, 9, 9, 10, 46, 10, 47, 46, 0, 11, 10, 10, 11, 47, 11, 48, 47, 0, 12, 11, 11, 12, 48, 12, 49, 48, 0, 13, 12, 12, 13, 49, 13, 50, 49, 0, 14, 13, 13, 14, 50, 14, 51, 50, 0, 15, 14, 14, 15, 51, 15, 52, 51, 0, 16, 15, 15, 16, 52, 16, 53, 52, 0, 17, 16, 16, 17, 53, 17, 54, 53, 0, 18, 17, 17, 18, 54, 18, 55, 54, 0, 19, 18, 18, 19, 55, 19, 56, 55, 0, 20, 19, 19, 20, 56, 20, 57, 56, 0, 21, 20, 20, 21, 57, 21, 58, 57, 0, 22, 21, 21, 22, 58, 22, 59, 58, 0, 23, 22, 22, 23, 59, 23, 60, 59, 0, 24, 23, 23, 24, 60, 24, 61, 60, 0, 25, 24, 24, 25, 61, 25, 62, 61, 0, 26, 25, 25, 26, 62, 26, 63, 62, 0, 27, 26, 26, 27, 63, 27, 64, 63, 0, 28, 27, 27, 28, 64, 28, 65, 64, 0, 29, 28, 28, 29, 65, 29, 66, 65, 0, 30, 29, 29, 30, 66, 30, 67, 66, 0, 31, 30, 30, 31, 67, 31, 68, 67, 0, 32, 31, 31, 32, 68, 32, 69, 68, 0, 33, 32, 32, 33, 69, 33, 70, 69, 0, 34, 33, 33, 34, 70, 34, 71, 70, 0, 35, 34, 34, 35, 71, 35, 72, 71, 0, 36, 35, 35, 36, 72, 36, 73, 72, 0, 1, 36, 36, 1, 73, 1, 38, 73, 37, 38, 39, 37, 39, 40, 37, 40, 41, 37, 41, 42, 37, 42, 43, 37, 43, 44, 37, 44, 45, 37, 45, 46, 37, 46, 47, 37, 47, 48, 37, 48, 49, 37, 49, 50, 37, 50, 51, 37, 51, 52, 37, 52, 53, 37, 53, 54, 37, 54, 55, 37, 55, 56, 37, 56, 57, 37, 57, 58, 37, 58, 59, 37, 59, 60, 37, 60, 61, 37, 61, 62, 37, 62, 63, 37, 63, 64, 37, 64, 65, 37, 65, 66, 37, 66, 67, 37, 67, 68, 37, 68, 69, 37, 69, 70, 37, 70, 71, 37, 71, 72, 37, 72, 73, 37, 73, 38], [ 0, 2, 1, 1, 2, 38, 2, 39, 38, 0, 3, 2, 2, 3, 39, 3, 40, 39, 0, 4, 3, 3, 4, 40, 4, 41, 40, 0, 5, 4, 4, 5, 41, 5, 42, 41, 0, 6, 5, 5, 6, 42, 6, 43, 42, 0, 7, 6, 6, 7, 43, 7, 44, 43, 0, 8, 7, 7, 8, 44, 8, 45, 44, 0, 9, 8, 8, 9, 45, 9, 46, 45, 0, 10, 9, 9, 10, 46, 10, 47, 46, 0, 11, 10, 10, 11, 47, 11, 48, 47, 0, 12, 11, 11, 12, 48, 12, 49, 48, 0, 13, 12, 12, 13, 49, 13, 50, 49, 0, 14, 13, 13, 14, 50, 14, 51, 50, 0, 15, 14, 14, 15, 51, 15, 52, 51, 0, 16, 15, 15, 16, 52, 16, 53, 52, 0, 17, 16, 16, 17, 53, 17, 54, 53, 0, 18, 17, 17, 18, 54, 18, 55, 54, 0, 19, 18, 18, 19, 55, 19, 56, 55, 0, 20, 19, 19, 20, 56, 20, 57, 56, 0, 21, 20, 20, 21, 57, 21, 58, 57, 0, 22, 21, 21, 22, 58, 22, 59, 58, 0, 23, 22, 22, 23, 59, 23, 60, 59, 0, 24, 23, 23, 24, 60, 24, 61, 60, 0, 25, 24, 24, 25, 61, 25, 62, 61, 0, 26, 25, 25, 26, 62, 26, 63, 62, 0, 27, 26, 26, 27, 63, 27, 64, 63, 0, 28, 27, 27, 28, 64, 28, 65, 64, 0, 29, 28, 28, 29, 65, 29, 66, 65, 0, 30, 29, 29, 30, 66, 30, 67, 66, 0, 31, 30, 30, 31, 67, 31, 68, 67, 0, 32, 31, 31, 32, 68, 32, 69, 68, 0, 33, 32, 32, 33, 69, 33, 70, 69, 0, 34, 33, 33, 34, 70, 34, 71, 70, 0, 35, 34, 34, 35, 71, 35, 72, 71, 0, 36, 35, 35, 36, 72, 36, 73, 72, 0, 1, 36, 36, 1, 73, 1, 38, 73, 37, 38, 39, 37, 39, 40, 37, 40, 41, 37, 41, 42, 37, 42, 43, 37, 43, 44, 37, 44, 45, 37, 45, 46, 37, 46, 47, 37, 47, 48, 37, 48, 49, 37, 49, 50, 37, 50, 51, 37, 51, 52, 37, 52, 53, 37, 53, 54, 37, 54, 55, 37, 55, 56, 37, 56, 57, 37, 57, 58, 37, 58, 59, 37, 59, 60, 37, 60, 61, 37, 61, 62, 37, 62, 63, 37, 63, 64, 37, 64, 65, 37, 65, 66, 37, 66, 67, 37, 67, 68, 37, 68, 69, 37, 69, 70, 37, 70, 71, 37, 71, 72, 37, 72, 73, 37, 73, 38], [ 0, 2, 1, 1, 2, 38, 2, 39, 38, 0, 3, 2, 2, 3, 39, 3, 40, 39, 0, 4, 3, 3, 4, 40, 4, 41, 40, 0, 5, 4, 4, 5, 41, 5, 42, 41, 0, 6, 5, 5, 6, 42, 6, 43, 42, 0, 7, 6, 6, 7, 43, 7, 44, 43, 0, 8, 7, 7, 8, 44, 8, 45, 44, 0, 9, 8, 8, 9, 45, 9, 46, 45, 0, 10, 9, 9, 10, 46, 10, 47, 46, 0, 11, 10, 10, 11, 47, 11, 48, 47, 0, 12, 11, 11, 12, 48, 12, 49, 48, 0, 13, 12, 12, 13, 49, 13, 50, 49, 0, 14, 13, 13, 14, 50, 14, 51, 50, 0, 15, 14, 14, 15, 51, 15, 52, 51, 0, 16, 15, 15, 16, 52, 16, 53, 52, 0, 17, 16, 16, 17, 53, 17, 54, 53, 0, 18, 17, 17, 18, 54, 18, 55, 54, 0, 19, 18, 18, 19, 55, 19, 56, 55, 0, 20, 19, 19, 20, 56, 20, 57, 56, 0, 21, 20, 20, 21, 57, 21, 58, 57, 0, 22, 21, 21, 22, 58, 22, 59, 58, 0, 23, 22, 22, 23, 59, 23, 60, 59, 0, 24, 23, 23, 24, 60, 24, 61, 60, 0, 25, 24, 24, 25, 61, 25, 62, 61, 0, 26, 25, 25, 26, 62, 26, 63, 62, 0, 27, 26, 26, 27, 63, 27, 64, 63, 0, 28, 27, 27, 28, 64, 28, 65, 64, 0, 29, 28, 28, 29, 65, 29, 66, 65, 0, 30, 29, 29, 30, 66, 30, 67, 66, 0, 31, 30, 30, 31, 67, 31, 68, 67, 0, 32, 31, 31, 32, 68, 32, 69, 68, 0, 33, 32, 32, 33, 69, 33, 70, 69, 0, 34, 33, 33, 34, 70, 34, 71, 70, 0, 35, 34, 34, 35, 71, 35, 72, 71, 0, 36, 35, 35, 36, 72, 36, 73, 72, 0, 1, 36, 36, 1, 73, 1, 38, 73, 37, 38, 39, 37, 39, 40, 37, 40, 41, 37, 41, 42, 37, 42, 43, 37, 43, 44, 37, 44, 45, 37, 45, 46, 37, 46, 47, 37, 47, 48, 37, 48, 49, 37, 49, 50, 37, 50, 51, 37, 51, 52, 37, 52, 53, 37, 53, 54, 37, 54, 55, 37, 55, 56, 37, 56, 57, 37, 57, 58, 37, 58, 59, 37, 59, 60, 37, 60, 61, 37, 61, 62, 37, 62, 63, 37, 63, 64, 37, 64, 65, 37, 65, 66, 37, 66, 67, 37, 67, 68, 37, 68, 69, 37, 69, 70, 37, 70, 71, 37, 71, 72, 37, 72, 73, 37, 73, 38], [ 0, 2, 1, 1, 2, 38, 2, 39, 38, 0, 3, 2, 2, 3, 39, 3, 40, 39, 0, 4, 3, 3, 4, 40, 4, 41, 40, 0, 5, 4, 4, 5, 41, 5, 42, 41, 0, 6, 5, 5, 6, 42, 6, 43, 42, 0, 7, 6, 6, 7, 43, 7, 44, 43, 0, 8, 7, 7, 8, 44, 8, 45, 44, 0, 9, 8, 8, 9, 45, 9, 46, 45, 0, 10, 9, 9, 10, 46, 10, 47, 46, 0, 11, 10, 10, 11, 47, 11, 48, 47, 0, 12, 11, 11, 12, 48, 12, 49, 48, 0, 13, 12, 12, 13, 49, 13, 50, 49, 0, 14, 13, 13, 14, 50, 14, 51, 50, 0, 15, 14, 14, 15, 51, 15, 52, 51, 0, 16, 15, 15, 16, 52, 16, 53, 52, 0, 17, 16, 16, 17, 53, 17, 54, 53, 0, 18, 17, 17, 18, 54, 18, 55, 54, 0, 19, 18, 18, 19, 55, 19, 56, 55, 0, 20, 19, 19, 20, 56, 20, 57, 56, 0, 21, 20, 20, 21, 57, 21, 58, 57, 0, 22, 21, 21, 22, 58, 22, 59, 58, 0, 23, 22, 22, 23, 59, 23, 60, 59, 0, 24, 23, 23, 24, 60, 24, 61, 60, 0, 25, 24, 24, 25, 61, 25, 62, 61, 0, 26, 25, 25, 26, 62, 26, 63, 62, 0, 27, 26, 26, 27, 63, 27, 64, 63, 0, 28, 27, 27, 28, 64, 28, 65, 64, 0, 29, 28, 28, 29, 65, 29, 66, 65, 0, 30, 29, 29, 30, 66, 30, 67, 66, 0, 31, 30, 30, 31, 67, 31, 68, 67, 0, 32, 31, 31, 32, 68, 32, 69, 68, 0, 33, 32, 32, 33, 69, 33, 70, 69, 0, 34, 33, 33, 34, 70, 34, 71, 70, 0, 35, 34, 34, 35, 71, 35, 72, 71, 0, 36, 35, 35, 36, 72, 36, 73, 72, 0, 1, 36, 36, 1, 73, 1, 38, 73, 37, 38, 39, 37, 39, 40, 37, 40, 41, 37, 41, 42, 37, 42, 43, 37, 43, 44, 37, 44, 45, 37, 45, 46, 37, 46, 47, 37, 47, 48, 37, 48, 49, 37, 49, 50, 37, 50, 51, 37, 51, 52, 37, 52, 53, 37, 53, 54, 37, 54, 55, 37, 55, 56, 37, 56, 57, 37, 57, 58, 37, 58, 59, 37, 59, 60, 37, 60, 61, 37, 61, 62, 37, 62, 63, 37, 63, 64, 37, 64, 65, 37, 65, 66, 37, 66, 67, 37, 67, 68, 37, 68, 69, 37, 69, 70, 37, 70, 71, 37, 71, 72, 37, 72, 73, 37, 73, 38], [ 0, 2, 1, 1, 2, 38, 2, 39, 38, 0, 3, 2, 2, 3, 39, 3, 40, 39, 0, 4, 3, 3, 4, 40, 4, 41, 40, 0, 5, 4, 4, 5, 41, 5, 42, 41, 0, 6, 5, 5, 6, 42, 6, 43, 42, 0, 7, 6, 6, 7, 43, 7, 44, 43, 0, 8, 7, 7, 8, 44, 8, 45, 44, 0, 9, 8, 8, 9, 45, 9, 46, 45, 0, 10, 9, 9, 10, 46, 10, 47, 46, 0, 11, 10, 10, 11, 47, 11, 48, 47, 0, 12, 11, 11, 12, 48, 12, 49, 48, 0, 13, 12, 12, 13, 49, 13, 50, 49, 0, 14, 13, 13, 14, 50, 14, 51, 50, 0, 15, 14, 14, 15, 51, 15, 52, 51, 0, 16, 15, 15, 16, 52, 16, 53, 52, 0, 17, 16, 16, 17, 53, 17, 54, 53, 0, 18, 17, 17, 18, 54, 18, 55, 54, 0, 19, 18, 18, 19, 55, 19, 56, 55, 0, 20, 19, 19, 20, 56, 20, 57, 56, 0, 21, 20, 20, 21, 57, 21, 58, 57, 0, 22, 21, 21, 22, 58, 22, 59, 58, 0, 23, 22, 22, 23, 59, 23, 60, 59, 0, 24, 23, 23, 24, 60, 24, 61, 60, 0, 25, 24, 24, 25, 61, 25, 62, 61, 0, 26, 25, 25, 26, 62, 26, 63, 62, 0, 27, 26, 26, 27, 63, 27, 64, 63, 0, 28, 27, 27, 28, 64, 28, 65, 64, 0, 29, 28, 28, 29, 65, 29, 66, 65, 0, 30, 29, 29, 30, 66, 30, 67, 66, 0, 31, 30, 30, 31, 67, 31, 68, 67, 0, 32, 31, 31, 32, 68, 32, 69, 68, 0, 33, 32, 32, 33, 69, 33, 70, 69, 0, 34, 33, 33, 34, 70, 34, 71, 70, 0, 35, 34, 34, 35, 71, 35, 72, 71, 0, 36, 35, 35, 36, 72, 36, 73, 72, 0, 1, 36, 36, 1, 73, 1, 38, 73, 37, 38, 39, 37, 39, 40, 37, 40, 41, 37, 41, 42, 37, 42, 43, 37, 43, 44, 37, 44, 45, 37, 45, 46, 37, 46, 47, 37, 47, 48, 37, 48, 49, 37, 49, 50, 37, 50, 51, 37, 51, 52, 37, 52, 53, 37, 53, 54, 37, 54, 55, 37, 55, 56, 37, 56, 57, 37, 57, 58, 37, 58, 59, 37, 59, 60, 37, 60, 61, 37, 61, 62, 37, 62, 63, 37, 63, 64, 37, 64, 65, 37, 65, 66, 37, 66, 67, 37, 67, 68, 37, 68, 69, 37, 69, 70, 37, 70, 71, 37, 71, 72, 37, 72, 73, 37, 73, 38], [ 0, 2, 1, 1, 2, 38, 2, 39, 38, 0, 3, 2, 2, 3, 39, 3, 40, 39, 0, 4, 3, 3, 4, 40, 4, 41, 40, 0, 5, 4, 4, 5, 41, 5, 42, 41, 0, 6, 5, 5, 6, 42, 6, 43, 42, 0, 7, 6, 6, 7, 43, 7, 44, 43, 0, 8, 7, 7, 8, 44, 8, 45, 44, 0, 9, 8, 8, 9, 45, 9, 46, 45, 0, 10, 9, 9, 10, 46, 10, 47, 46, 0, 11, 10, 10, 11, 47, 11, 48, 47, 0, 12, 11, 11, 12, 48, 12, 49, 48, 0, 13, 12, 12, 13, 49, 13, 50, 49, 0, 14, 13, 13, 14, 50, 14, 51, 50, 0, 15, 14, 14, 15, 51, 15, 52, 51, 0, 16, 15, 15, 16, 52, 16, 53, 52, 0, 17, 16, 16, 17, 53, 17, 54, 53, 0, 18, 17, 17, 18, 54, 18, 55, 54, 0, 19, 18, 18, 19, 55, 19, 56, 55, 0, 20, 19, 19, 20, 56, 20, 57, 56, 0, 21, 20, 20, 21, 57, 21, 58, 57, 0, 22, 21, 21, 22, 58, 22, 59, 58, 0, 23, 22, 22, 23, 59, 23, 60, 59, 0, 24, 23, 23, 24, 60, 24, 61, 60, 0, 25, 24, 24, 25, 61, 25, 62, 61, 0, 26, 25, 25, 26, 62, 26, 63, 62, 0, 27, 26, 26, 27, 63, 27, 64, 63, 0, 28, 27, 27, 28, 64, 28, 65, 64, 0, 29, 28, 28, 29, 65, 29, 66, 65, 0, 30, 29, 29, 30, 66, 30, 67, 66, 0, 31, 30, 30, 31, 67, 31, 68, 67, 0, 32, 31, 31, 32, 68, 32, 69, 68, 0, 33, 32, 32, 33, 69, 33, 70, 69, 0, 34, 33, 33, 34, 70, 34, 71, 70, 0, 35, 34, 34, 35, 71, 35, 72, 71, 0, 36, 35, 35, 36, 72, 36, 73, 72, 0, 1, 36, 36, 1, 73, 1, 38, 73, 37, 38, 39, 37, 39, 40, 37, 40, 41, 37, 41, 42, 37, 42, 43, 37, 43, 44, 37, 44, 45, 37, 45, 46, 37, 46, 47, 37, 47, 48, 37, 48, 49, 37, 49, 50, 37, 50, 51, 37, 51, 52, 37, 52, 53, 37, 53, 54, 37, 54, 55, 37, 55, 56, 37, 56, 57, 37, 57, 58, 37, 58, 59, 37, 59, 60, 37, 60, 61, 37, 61, 62, 37, 62, 63, 37, 63, 64, 37, 64, 65, 37, 65, 66, 37, 66, 67, 37, 67, 68, 37, 68, 69, 37, 69, 70, 37, 70, 71, 37, 71, 72, 37, 72, 73, 37, 73, 38], [ 0, 2, 1, 1, 2, 38, 2, 39, 38, 0, 3, 2, 2, 3, 39, 3, 40, 39, 0, 4, 3, 3, 4, 40, 4, 41, 40, 0, 5, 4, 4, 5, 41, 5, 42, 41, 0, 6, 5, 5, 6, 42, 6, 43, 42, 0, 7, 6, 6, 7, 43, 7, 44, 43, 0, 8, 7, 7, 8, 44, 8, 45, 44, 0, 9, 8, 8, 9, 45, 9, 46, 45, 0, 10, 9, 9, 10, 46, 10, 47, 46, 0, 11, 10, 10, 11, 47, 11, 48, 47, 0, 12, 11, 11, 12, 48, 12, 49, 48, 0, 13, 12, 12, 13, 49, 13, 50, 49, 0, 14, 13, 13, 14, 50, 14, 51, 50, 0, 15, 14, 14, 15, 51, 15, 52, 51, 0, 16, 15, 15, 16, 52, 16, 53, 52, 0, 17, 16, 16, 17, 53, 17, 54, 53, 0, 18, 17, 17, 18, 54, 18, 55, 54, 0, 19, 18, 18, 19, 55, 19, 56, 55, 0, 20, 19, 19, 20, 56, 20, 57, 56, 0, 21, 20, 20, 21, 57, 21, 58, 57, 0, 22, 21, 21, 22, 58, 22, 59, 58, 0, 23, 22, 22, 23, 59, 23, 60, 59, 0, 24, 23, 23, 24, 60, 24, 61, 60, 0, 25, 24, 24, 25, 61, 25, 62, 61, 0, 26, 25, 25, 26, 62, 26, 63, 62, 0, 27, 26, 26, 27, 63, 27, 64, 63, 0, 28, 27, 27, 28, 64, 28, 65, 64, 0, 29, 28, 28, 29, 65, 29, 66, 65, 0, 30, 29, 29, 30, 66, 30, 67, 66, 0, 31, 30, 30, 31, 67, 31, 68, 67, 0, 32, 31, 31, 32, 68, 32, 69, 68, 0, 33, 32, 32, 33, 69, 33, 70, 69, 0, 34, 33, 33, 34, 70, 34, 71, 70, 0, 35, 34, 34, 35, 71, 35, 72, 71, 0, 36, 35, 35, 36, 72, 36, 73, 72, 0, 1, 36, 36, 1, 73, 1, 38, 73, 37, 38, 39, 37, 39, 40, 37, 40, 41, 37, 41, 42, 37, 42, 43, 37, 43, 44, 37, 44, 45, 37, 45, 46, 37, 46, 47, 37, 47, 48, 37, 48, 49, 37, 49, 50, 37, 50, 51, 37, 51, 52, 37, 52, 53, 37, 53, 54, 37, 54, 55, 37, 55, 56, 37, 56, 57, 37, 57, 58, 37, 58, 59, 37, 59, 60, 37, 60, 61, 37, 61, 62, 37, 62, 63, 37, 63, 64, 37, 64, 65, 37, 65, 66, 37, 66, 67, 37, 67, 68, 37, 68, 69, 37, 69, 70, 37, 70, 71, 37, 71, 72, 37, 72, 73, 37, 73, 38], [ 2, 1, 0, 2, 4, 3, 2, 5, 4, 2, 0, 37, 2, 37, 5, 37, 36, 5, 5, 36, 6, 21, 20, 19, 21, 23, 22, 21, 24, 23, 21, 19, 18, 21, 18, 24, 18, 17, 24, 24, 17, 25, 40, 38, 39, 40, 41, 42, 40, 42, 43, 40, 75, 38, 40, 43, 75, 75, 43, 74, 43, 44, 74, 59, 57, 58, 59, 60, 61, 59, 61, 62, 59, 56, 57, 59, 62, 56, 56, 62, 55, 62, 63, 55, 6, 36, 7, 36, 35, 7, 44, 45, 74, 74, 45, 73, 7, 35, 8, 35, 34, 8, 45, 46, 73, 73, 46, 72, 8, 34, 9, 34, 33, 9, 46, 47, 72, 72, 47, 71, 9, 33, 10, 33, 32, 10, 47, 48, 71, 71, 48, 70, 10, 32, 11, 32, 31, 11, 48, 49, 70, 70, 49, 69, 11, 31, 12, 31, 30, 12, 49, 50, 69, 69, 50, 68, 12, 30, 13, 30, 29, 13, 50, 51, 68, 68, 51, 67, 13, 29, 14, 29, 28, 14, 51, 52, 67, 67, 52, 66, 14, 28, 15, 28, 27, 15, 52, 53, 66, 66, 53, 65, 15, 27, 16, 27, 26, 16, 53, 54, 65, 65, 54, 64, 16, 26, 17, 26, 25, 17, 54, 55, 64, 64, 55, 63, 0, 1, 39, 0, 39, 38, 1, 2, 40, 1, 40, 39, 2, 3, 41, 2, 41, 40, 3, 4, 42, 3, 42, 41, 6, 7, 45, 6, 45, 44, 7, 8, 46, 7, 46, 45, 8, 9, 47, 8, 47, 46, 9, 10, 48, 9, 48, 47, 10, 11, 49, 10, 49, 48, 11, 12, 50, 11, 50, 49, 12, 13, 51, 12, 51, 50, 13, 14, 52, 13, 52, 51, 14, 15, 53, 14, 53, 52, 15, 16, 54, 15, 54, 53, 16, 17, 55, 16, 55, 54, 19, 20, 58, 19, 58, 57, 20, 21, 59, 20, 59, 58, 21, 22, 60, 21, 60, 59, 22, 23, 61, 22, 61, 60, 25, 26, 64, 25, 64, 63, 26, 27, 65, 26, 65, 64, 27, 28, 66, 27, 66, 65, 28, 29, 67, 28, 67, 66, 29, 30, 68, 29, 68, 67, 30, 31, 69, 30, 69, 68, 31, 32, 70, 31, 70, 69, 32, 33, 71, 32, 71, 70, 33, 34, 72, 33, 72, 71, 34, 35, 73, 34, 73, 72, 35, 36, 74, 35, 74, 73, 37, 0, 36, 76, 36, 0, 76, 74, 36, 76, 38, 74, 75, 74, 38, 4, 5, 6, 4, 6, 77, 77, 6, 44, 42, 77, 44, 42, 44, 43, 18, 19, 17, 17, 19, 78, 17, 78, 55, 55, 78, 57, 56, 55, 57, 25, 23, 24, 25, 79, 23, 79, 25, 63, 61, 79, 63, 61, 63, 62], [ 0, 2, 1, 1, 2, 38, 2, 39, 38, 0, 3, 2, 2, 3, 39, 3, 40, 39, 0, 4, 3, 3, 4, 40, 4, 41, 40, 0, 5, 4, 4, 5, 41, 5, 42, 41, 0, 6, 5, 5, 6, 42, 6, 43, 42, 0, 7, 6, 6, 7, 43, 7, 44, 43, 0, 8, 7, 7, 8, 44, 8, 45, 44, 0, 9, 8, 8, 9, 45, 9, 46, 45, 0, 10, 9, 9, 10, 46, 10, 47, 46, 0, 11, 10, 10, 11, 47, 11, 48, 47, 0, 12, 11, 11, 12, 48, 12, 49, 48, 0, 13, 12, 12, 13, 49, 13, 50, 49, 0, 14, 13, 13, 14, 50, 14, 51, 50, 0, 15, 14, 14, 15, 51, 15, 52, 51, 0, 16, 15, 15, 16, 52, 16, 53, 52, 0, 17, 16, 16, 17, 53, 17, 54, 53, 0, 18, 17, 17, 18, 54, 18, 55, 54, 0, 19, 18, 18, 19, 55, 19, 56, 55, 0, 20, 19, 19, 20, 56, 20, 57, 56, 0, 21, 20, 20, 21, 57, 21, 58, 57, 0, 22, 21, 21, 22, 58, 22, 59, 58, 0, 23, 22, 22, 23, 59, 23, 60, 59, 0, 24, 23, 23, 24, 60, 24, 61, 60, 0, 25, 24, 24, 25, 61, 25, 62, 61, 0, 26, 25, 25, 26, 62, 26, 63, 62, 0, 27, 26, 26, 27, 63, 27, 64, 63, 0, 28, 27, 27, 28, 64, 28, 65, 64, 0, 29, 28, 28, 29, 65, 29, 66, 65, 0, 30, 29, 29, 30, 66, 30, 67, 66, 0, 31, 30, 30, 31, 67, 31, 68, 67, 0, 32, 31, 31, 32, 68, 32, 69, 68, 0, 33, 32, 32, 33, 69, 33, 70, 69, 0, 34, 33, 33, 34, 70, 34, 71, 70, 0, 35, 34, 34, 35, 71, 35, 72, 71, 0, 36, 35, 35, 36, 72, 36, 73, 72, 0, 1, 36, 36, 1, 73, 1, 38, 73, 37, 38, 39, 37, 39, 40, 37, 40, 41, 37, 41, 42, 37, 42, 43, 37, 43, 44, 37, 44, 45, 37, 45, 46, 37, 46, 47, 37, 47, 48, 37, 48, 49, 37, 49, 50, 37, 50, 51, 37, 51, 52, 37, 52, 53, 37, 53, 54, 37, 54, 55, 37, 55, 56, 37, 56, 57, 37, 57, 58, 37, 58, 59, 37, 59, 60, 37, 60, 61, 37, 61, 62, 37, 62, 63, 37, 63, 64, 37, 64, 65, 37, 65, 66, 37, 66, 67, 37, 67, 68, 37, 68, 69, 37, 69, 70, 37, 70, 71, 37, 71, 72, 37, 72, 73, 37, 73, 38], [ 0, 2, 1, 1, 2, 38, 2, 39, 38, 0, 3, 2, 2, 3, 39, 3, 40, 39, 0, 4, 3, 3, 4, 40, 4, 41, 40, 0, 5, 4, 4, 5, 41, 5, 42, 41, 0, 6, 5, 5, 6, 42, 6, 43, 42, 0, 7, 6, 6, 7, 43, 7, 44, 43, 0, 8, 7, 7, 8, 44, 8, 45, 44, 0, 9, 8, 8, 9, 45, 9, 46, 45, 0, 10, 9, 9, 10, 46, 10, 47, 46, 0, 11, 10, 10, 11, 47, 11, 48, 47, 0, 12, 11, 11, 12, 48, 12, 49, 48, 0, 13, 12, 12, 13, 49, 13, 50, 49, 0, 14, 13, 13, 14, 50, 14, 51, 50, 0, 15, 14, 14, 15, 51, 15, 52, 51, 0, 16, 15, 15, 16, 52, 16, 53, 52, 0, 17, 16, 16, 17, 53, 17, 54, 53, 0, 18, 17, 17, 18, 54, 18, 55, 54, 0, 19, 18, 18, 19, 55, 19, 56, 55, 0, 20, 19, 19, 20, 56, 20, 57, 56, 0, 21, 20, 20, 21, 57, 21, 58, 57, 0, 22, 21, 21, 22, 58, 22, 59, 58, 0, 23, 22, 22, 23, 59, 23, 60, 59, 0, 24, 23, 23, 24, 60, 24, 61, 60, 0, 25, 24, 24, 25, 61, 25, 62, 61, 0, 26, 25, 25, 26, 62, 26, 63, 62, 0, 27, 26, 26, 27, 63, 27, 64, 63, 0, 28, 27, 27, 28, 64, 28, 65, 64, 0, 29, 28, 28, 29, 65, 29, 66, 65, 0, 30, 29, 29, 30, 66, 30, 67, 66, 0, 31, 30, 30, 31, 67, 31, 68, 67, 0, 32, 31, 31, 32, 68, 32, 69, 68, 0, 33, 32, 32, 33, 69, 33, 70, 69, 0, 34, 33, 33, 34, 70, 34, 71, 70, 0, 35, 34, 34, 35, 71, 35, 72, 71, 0, 36, 35, 35, 36, 72, 36, 73, 72, 0, 1, 36, 36, 1, 73, 1, 38, 73, 37, 38, 39, 37, 39, 40, 37, 40, 41, 37, 41, 42, 37, 42, 43, 37, 43, 44, 37, 44, 45, 37, 45, 46, 37, 46, 47, 37, 47, 48, 37, 48, 49, 37, 49, 50, 37, 50, 51, 37, 51, 52, 37, 52, 53, 37, 53, 54, 37, 54, 55, 37, 55, 56, 37, 56, 57, 37, 57, 58, 37, 58, 59, 37, 59, 60, 37, 60, 61, 37, 61, 62, 37, 62, 63, 37, 63, 64, 37, 64, 65, 37, 65, 66, 37, 66, 67, 37, 67, 68, 37, 68, 69, 37, 69, 70, 37, 70, 71, 37, 71, 72, 37, 72, 73, 37, 73, 38], [ 0, 2, 1, 1, 2, 38, 2, 39, 38, 0, 3, 2, 2, 3, 39, 3, 40, 39, 0, 4, 3, 3, 4, 40, 4, 41, 40, 0, 5, 4, 4, 5, 41, 5, 42, 41, 0, 6, 5, 5, 6, 42, 6, 43, 42, 0, 7, 6, 6, 7, 43, 7, 44, 43, 0, 8, 7, 7, 8, 44, 8, 45, 44, 0, 9, 8, 8, 9, 45, 9, 46, 45, 0, 10, 9, 9, 10, 46, 10, 47, 46, 0, 11, 10, 10, 11, 47, 11, 48, 47, 0, 12, 11, 11, 12, 48, 12, 49, 48, 0, 13, 12, 12, 13, 49, 13, 50, 49, 0, 14, 13, 13, 14, 50, 14, 51, 50, 0, 15, 14, 14, 15, 51, 15, 52, 51, 0, 16, 15, 15, 16, 52, 16, 53, 52, 0, 17, 16, 16, 17, 53, 17, 54, 53, 0, 18, 17, 17, 18, 54, 18, 55, 54, 0, 19, 18, 18, 19, 55, 19, 56, 55, 0, 20, 19, 19, 20, 56, 20, 57, 56, 0, 21, 20, 20, 21, 57, 21, 58, 57, 0, 22, 21, 21, 22, 58, 22, 59, 58, 0, 23, 22, 22, 23, 59, 23, 60, 59, 0, 24, 23, 23, 24, 60, 24, 61, 60, 0, 25, 24, 24, 25, 61, 25, 62, 61, 0, 26, 25, 25, 26, 62, 26, 63, 62, 0, 27, 26, 26, 27, 63, 27, 64, 63, 0, 28, 27, 27, 28, 64, 28, 65, 64, 0, 29, 28, 28, 29, 65, 29, 66, 65, 0, 30, 29, 29, 30, 66, 30, 67, 66, 0, 31, 30, 30, 31, 67, 31, 68, 67, 0, 32, 31, 31, 32, 68, 32, 69, 68, 0, 33, 32, 32, 33, 69, 33, 70, 69, 0, 34, 33, 33, 34, 70, 34, 71, 70, 0, 35, 34, 34, 35, 71, 35, 72, 71, 0, 36, 35, 35, 36, 72, 36, 73, 72, 0, 1, 36, 36, 1, 73, 1, 38, 73, 37, 38, 39, 37, 39, 40, 37, 40, 41, 37, 41, 42, 37, 42, 43, 37, 43, 44, 37, 44, 45, 37, 45, 46, 37, 46, 47, 37, 47, 48, 37, 48, 49, 37, 49, 50, 37, 50, 51, 37, 51, 52, 37, 52, 53, 37, 53, 54, 37, 54, 55, 37, 55, 56, 37, 56, 57, 37, 57, 58, 37, 58, 59, 37, 59, 60, 37, 60, 61, 37, 61, 62, 37, 62, 63, 37, 63, 64, 37, 64, 65, 37, 65, 66, 37, 66, 67, 37, 67, 68, 37, 68, 69, 37, 69, 70, 37, 70, 71, 37, 71, 72, 37, 72, 73, 37, 73, 38], [ 0, 2, 1, 1, 2, 38, 2, 39, 38, 0, 3, 2, 2, 3, 39, 3, 40, 39, 0, 4, 3, 3, 4, 40, 4, 41, 40, 0, 5, 4, 4, 5, 41, 5, 42, 41, 0, 6, 5, 5, 6, 42, 6, 43, 42, 0, 7, 6, 6, 7, 43, 7, 44, 43, 0, 8, 7, 7, 8, 44, 8, 45, 44, 0, 9, 8, 8, 9, 45, 9, 46, 45, 0, 10, 9, 9, 10, 46, 10, 47, 46, 0, 11, 10, 10, 11, 47, 11, 48, 47, 0, 12, 11, 11, 12, 48, 12, 49, 48, 0, 13, 12, 12, 13, 49, 13, 50, 49, 0, 14, 13, 13, 14, 50, 14, 51, 50, 0, 15, 14, 14, 15, 51, 15, 52, 51, 0, 16, 15, 15, 16, 52, 16, 53, 52, 0, 17, 16, 16, 17, 53, 17, 54, 53, 0, 18, 17, 17, 18, 54, 18, 55, 54, 0, 19, 18, 18, 19, 55, 19, 56, 55, 0, 20, 19, 19, 20, 56, 20, 57, 56, 0, 21, 20, 20, 21, 57, 21, 58, 57, 0, 22, 21, 21, 22, 58, 22, 59, 58, 0, 23, 22, 22, 23, 59, 23, 60, 59, 0, 24, 23, 23, 24, 60, 24, 61, 60, 0, 25, 24, 24, 25, 61, 25, 62, 61, 0, 26, 25, 25, 26, 62, 26, 63, 62, 0, 27, 26, 26, 27, 63, 27, 64, 63, 0, 28, 27, 27, 28, 64, 28, 65, 64, 0, 29, 28, 28, 29, 65, 29, 66, 65, 0, 30, 29, 29, 30, 66, 30, 67, 66, 0, 31, 30, 30, 31, 67, 31, 68, 67, 0, 32, 31, 31, 32, 68, 32, 69, 68, 0, 33, 32, 32, 33, 69, 33, 70, 69, 0, 34, 33, 33, 34, 70, 34, 71, 70, 0, 35, 34, 34, 35, 71, 35, 72, 71, 0, 36, 35, 35, 36, 72, 36, 73, 72, 0, 1, 36, 36, 1, 73, 1, 38, 73, 37, 38, 39, 37, 39, 40, 37, 40, 41, 37, 41, 42, 37, 42, 43, 37, 43, 44, 37, 44, 45, 37, 45, 46, 37, 46, 47, 37, 47, 48, 37, 48, 49, 37, 49, 50, 37, 50, 51, 37, 51, 52, 37, 52, 53, 37, 53, 54, 37, 54, 55, 37, 55, 56, 37, 56, 57, 37, 57, 58, 37, 58, 59, 37, 59, 60, 37, 60, 61, 37, 61, 62, 37, 62, 63, 37, 63, 64, 37, 64, 65, 37, 65, 66, 37, 66, 67, 37, 67, 68, 37, 68, 69, 37, 69, 70, 37, 70, 71, 37, 71, 72, 37, 72, 73, 37, 73, 38], [ 2, 1, 0, 2, 4, 3, 2, 5, 4, 2, 0, 37, 2, 37, 5, 37, 36, 5, 5, 36, 6, 21, 20, 19, 21, 23, 22, 21, 24, 23, 21, 19, 18, 21, 18, 24, 18, 17, 24, 24, 17, 25, 40, 38, 39, 40, 41, 42, 40, 42, 43, 40, 75, 38, 40, 43, 75, 75, 43, 74, 43, 44, 74, 59, 57, 58, 59, 60, 61, 59, 61, 62, 59, 56, 57, 59, 62, 56, 56, 62, 55, 62, 63, 55, 6, 36, 7, 36, 35, 7, 44, 45, 74, 74, 45, 73, 7, 35, 8, 35, 34, 8, 45, 46, 73, 73, 46, 72, 8, 34, 9, 34, 33, 9, 46, 47, 72, 72, 47, 71, 9, 33, 10, 33, 32, 10, 47, 48, 71, 71, 48, 70, 10, 32, 11, 32, 31, 11, 48, 49, 70, 70, 49, 69, 11, 31, 12, 31, 30, 12, 49, 50, 69, 69, 50, 68, 12, 30, 13, 30, 29, 13, 50, 51, 68, 68, 51, 67, 13, 29, 14, 29, 28, 14, 51, 52, 67, 67, 52, 66, 14, 28, 15, 28, 27, 15, 52, 53, 66, 66, 53, 65, 15, 27, 16, 27, 26, 16, 53, 54, 65, 65, 54, 64, 16, 26, 17, 26, 25, 17, 54, 55, 64, 64, 55, 63, 0, 1, 39, 0, 39, 38, 1, 2, 40, 1, 40, 39, 2, 3, 41, 2, 41, 40, 3, 4, 42, 3, 42, 41, 6, 7, 45, 6, 45, 44, 7, 8, 46, 7, 46, 45, 8, 9, 47, 8, 47, 46, 9, 10, 48, 9, 48, 47, 10, 11, 49, 10, 49, 48, 11, 12, 50, 11, 50, 49, 12, 13, 51, 12, 51, 50, 13, 14, 52, 13, 52, 51, 14, 15, 53, 14, 53, 52, 15, 16, 54, 15, 54, 53, 16, 17, 55, 16, 55, 54, 19, 20, 58, 19, 58, 57, 20, 21, 59, 20, 59, 58, 21, 22, 60, 21, 60, 59, 22, 23, 61, 22, 61, 60, 25, 26, 64, 25, 64, 63, 26, 27, 65, 26, 65, 64, 27, 28, 66, 27, 66, 65, 28, 29, 67, 28, 67, 66, 29, 30, 68, 29, 68, 67, 30, 31, 69, 30, 69, 68, 31, 32, 70, 31, 70, 69, 32, 33, 71, 32, 71, 70, 33, 34, 72, 33, 72, 71, 34, 35, 73, 34, 73, 72, 35, 36, 74, 35, 74, 73, 37, 0, 36, 76, 36, 0, 76, 74, 36, 76, 38, 74, 75, 74, 38, 4, 5, 6, 4, 6, 77, 77, 6, 44, 42, 77, 44, 42, 44, 43, 18, 19, 17, 17, 19, 78, 17, 78, 55, 55, 78, 57, 56, 55, 57, 25, 23, 24, 25, 79, 23, 79, 25, 63, 61, 79, 63, 61, 63, 62], [ 0, 2, 1, 1, 2, 38, 2, 39, 38, 0, 3, 2, 2, 3, 39, 3, 40, 39, 0, 4, 3, 3, 4, 40, 4, 41, 40, 0, 5, 4, 4, 5, 41, 5, 42, 41, 0, 6, 5, 5, 6, 42, 6, 43, 42, 0, 7, 6, 6, 7, 43, 7, 44, 43, 0, 8, 7, 7, 8, 44, 8, 45, 44, 0, 9, 8, 8, 9, 45, 9, 46, 45, 0, 10, 9, 9, 10, 46, 10, 47, 46, 0, 11, 10, 10, 11, 47, 11, 48, 47, 0, 12, 11, 11, 12, 48, 12, 49, 48, 0, 13, 12, 12, 13, 49, 13, 50, 49, 0, 14, 13, 13, 14, 50, 14, 51, 50, 0, 15, 14, 14, 15, 51, 15, 52, 51, 0, 16, 15, 15, 16, 52, 16, 53, 52, 0, 17, 16, 16, 17, 53, 17, 54, 53, 0, 18, 17, 17, 18, 54, 18, 55, 54, 0, 19, 18, 18, 19, 55, 19, 56, 55, 0, 20, 19, 19, 20, 56, 20, 57, 56, 0, 21, 20, 20, 21, 57, 21, 58, 57, 0, 22, 21, 21, 22, 58, 22, 59, 58, 0, 23, 22, 22, 23, 59, 23, 60, 59, 0, 24, 23, 23, 24, 60, 24, 61, 60, 0, 25, 24, 24, 25, 61, 25, 62, 61, 0, 26, 25, 25, 26, 62, 26, 63, 62, 0, 27, 26, 26, 27, 63, 27, 64, 63, 0, 28, 27, 27, 28, 64, 28, 65, 64, 0, 29, 28, 28, 29, 65, 29, 66, 65, 0, 30, 29, 29, 30, 66, 30, 67, 66, 0, 31, 30, 30, 31, 67, 31, 68, 67, 0, 32, 31, 31, 32, 68, 32, 69, 68, 0, 33, 32, 32, 33, 69, 33, 70, 69, 0, 34, 33, 33, 34, 70, 34, 71, 70, 0, 35, 34, 34, 35, 71, 35, 72, 71, 0, 36, 35, 35, 36, 72, 36, 73, 72, 0, 1, 36, 36, 1, 73, 1, 38, 73, 37, 38, 39, 37, 39, 40, 37, 40, 41, 37, 41, 42, 37, 42, 43, 37, 43, 44, 37, 44, 45, 37, 45, 46, 37, 46, 47, 37, 47, 48, 37, 48, 49, 37, 49, 50, 37, 50, 51, 37, 51, 52, 37, 52, 53, 37, 53, 54, 37, 54, 55, 37, 55, 56, 37, 56, 57, 37, 57, 58, 37, 58, 59, 37, 59, 60, 37, 60, 61, 37, 61, 62, 37, 62, 63, 37, 63, 64, 37, 64, 65, 37, 65, 66, 37, 66, 67, 37, 67, 68, 37, 68, 69, 37, 69, 70, 37, 70, 71, 37, 71, 72, 37, 72, 73, 37, 73, 38], [ 0, 2, 1, 1, 2, 38, 2, 39, 38, 0, 3, 2, 2, 3, 39, 3, 40, 39, 0, 4, 3, 3, 4, 40, 4, 41, 40, 0, 5, 4, 4, 5, 41, 5, 42, 41, 0, 6, 5, 5, 6, 42, 6, 43, 42, 0, 7, 6, 6, 7, 43, 7, 44, 43, 0, 8, 7, 7, 8, 44, 8, 45, 44, 0, 9, 8, 8, 9, 45, 9, 46, 45, 0, 10, 9, 9, 10, 46, 10, 47, 46, 0, 11, 10, 10, 11, 47, 11, 48, 47, 0, 12, 11, 11, 12, 48, 12, 49, 48, 0, 13, 12, 12, 13, 49, 13, 50, 49, 0, 14, 13, 13, 14, 50, 14, 51, 50, 0, 15, 14, 14, 15, 51, 15, 52, 51, 0, 16, 15, 15, 16, 52, 16, 53, 52, 0, 17, 16, 16, 17, 53, 17, 54, 53, 0, 18, 17, 17, 18, 54, 18, 55, 54, 0, 19, 18, 18, 19, 55, 19, 56, 55, 0, 20, 19, 19, 20, 56, 20, 57, 56, 0, 21, 20, 20, 21, 57, 21, 58, 57, 0, 22, 21, 21, 22, 58, 22, 59, 58, 0, 23, 22, 22, 23, 59, 23, 60, 59, 0, 24, 23, 23, 24, 60, 24, 61, 60, 0, 25, 24, 24, 25, 61, 25, 62, 61, 0, 26, 25, 25, 26, 62, 26, 63, 62, 0, 27, 26, 26, 27, 63, 27, 64, 63, 0, 28, 27, 27, 28, 64, 28, 65, 64, 0, 29, 28, 28, 29, 65, 29, 66, 65, 0, 30, 29, 29, 30, 66, 30, 67, 66, 0, 31, 30, 30, 31, 67, 31, 68, 67, 0, 32, 31, 31, 32, 68, 32, 69, 68, 0, 33, 32, 32, 33, 69, 33, 70, 69, 0, 34, 33, 33, 34, 70, 34, 71, 70, 0, 35, 34, 34, 35, 71, 35, 72, 71, 0, 36, 35, 35, 36, 72, 36, 73, 72, 0, 1, 36, 36, 1, 73, 1, 38, 73, 37, 38, 39, 37, 39, 40, 37, 40, 41, 37, 41, 42, 37, 42, 43, 37, 43, 44, 37, 44, 45, 37, 45, 46, 37, 46, 47, 37, 47, 48, 37, 48, 49, 37, 49, 50, 37, 50, 51, 37, 51, 52, 37, 52, 53, 37, 53, 54, 37, 54, 55, 37, 55, 56, 37, 56, 57, 37, 57, 58, 37, 58, 59, 37, 59, 60, 37, 60, 61, 37, 61, 62, 37, 62, 63, 37, 63, 64, 37, 64, 65, 37, 65, 66, 37, 66, 67, 37, 67, 68, 37, 68, 69, 37, 69, 70, 37, 70, 71, 37, 71, 72, 37, 72, 73, 37, 73, 38], [ 0, 2, 1, 1, 2, 38, 2, 39, 38, 0, 3, 2, 2, 3, 39, 3, 40, 39, 0, 4, 3, 3, 4, 40, 4, 41, 40, 0, 5, 4, 4, 5, 41, 5, 42, 41, 0, 6, 5, 5, 6, 42, 6, 43, 42, 0, 7, 6, 6, 7, 43, 7, 44, 43, 0, 8, 7, 7, 8, 44, 8, 45, 44, 0, 9, 8, 8, 9, 45, 9, 46, 45, 0, 10, 9, 9, 10, 46, 10, 47, 46, 0, 11, 10, 10, 11, 47, 11, 48, 47, 0, 12, 11, 11, 12, 48, 12, 49, 48, 0, 13, 12, 12, 13, 49, 13, 50, 49, 0, 14, 13, 13, 14, 50, 14, 51, 50, 0, 15, 14, 14, 15, 51, 15, 52, 51, 0, 16, 15, 15, 16, 52, 16, 53, 52, 0, 17, 16, 16, 17, 53, 17, 54, 53, 0, 18, 17, 17, 18, 54, 18, 55, 54, 0, 19, 18, 18, 19, 55, 19, 56, 55, 0, 20, 19, 19, 20, 56, 20, 57, 56, 0, 21, 20, 20, 21, 57, 21, 58, 57, 0, 22, 21, 21, 22, 58, 22, 59, 58, 0, 23, 22, 22, 23, 59, 23, 60, 59, 0, 24, 23, 23, 24, 60, 24, 61, 60, 0, 25, 24, 24, 25, 61, 25, 62, 61, 0, 26, 25, 25, 26, 62, 26, 63, 62, 0, 27, 26, 26, 27, 63, 27, 64, 63, 0, 28, 27, 27, 28, 64, 28, 65, 64, 0, 29, 28, 28, 29, 65, 29, 66, 65, 0, 30, 29, 29, 30, 66, 30, 67, 66, 0, 31, 30, 30, 31, 67, 31, 68, 67, 0, 32, 31, 31, 32, 68, 32, 69, 68, 0, 33, 32, 32, 33, 69, 33, 70, 69, 0, 34, 33, 33, 34, 70, 34, 71, 70, 0, 35, 34, 34, 35, 71, 35, 72, 71, 0, 36, 35, 35, 36, 72, 36, 73, 72, 0, 1, 36, 36, 1, 73, 1, 38, 73, 37, 38, 39, 37, 39, 40, 37, 40, 41, 37, 41, 42, 37, 42, 43, 37, 43, 44, 37, 44, 45, 37, 45, 46, 37, 46, 47, 37, 47, 48, 37, 48, 49, 37, 49, 50, 37, 50, 51, 37, 51, 52, 37, 52, 53, 37, 53, 54, 37, 54, 55, 37, 55, 56, 37, 56, 57, 37, 57, 58, 37, 58, 59, 37, 59, 60, 37, 60, 61, 37, 61, 62, 37, 62, 63, 37, 63, 64, 37, 64, 65, 37, 65, 66, 37, 66, 67, 37, 67, 68, 37, 68, 69, 37, 69, 70, 37, 70, 71, 37, 71, 72, 37, 72, 73, 37, 73, 38], [ 0, 2, 1, 1, 2, 38, 2, 39, 38, 0, 3, 2, 2, 3, 39, 3, 40, 39, 0, 4, 3, 3, 4, 40, 4, 41, 40, 0, 5, 4, 4, 5, 41, 5, 42, 41, 0, 6, 5, 5, 6, 42, 6, 43, 42, 0, 7, 6, 6, 7, 43, 7, 44, 43, 0, 8, 7, 7, 8, 44, 8, 45, 44, 0, 9, 8, 8, 9, 45, 9, 46, 45, 0, 10, 9, 9, 10, 46, 10, 47, 46, 0, 11, 10, 10, 11, 47, 11, 48, 47, 0, 12, 11, 11, 12, 48, 12, 49, 48, 0, 13, 12, 12, 13, 49, 13, 50, 49, 0, 14, 13, 13, 14, 50, 14, 51, 50, 0, 15, 14, 14, 15, 51, 15, 52, 51, 0, 16, 15, 15, 16, 52, 16, 53, 52, 0, 17, 16, 16, 17, 53, 17, 54, 53, 0, 18, 17, 17, 18, 54, 18, 55, 54, 0, 19, 18, 18, 19, 55, 19, 56, 55, 0, 20, 19, 19, 20, 56, 20, 57, 56, 0, 21, 20, 20, 21, 57, 21, 58, 57, 0, 22, 21, 21, 22, 58, 22, 59, 58, 0, 23, 22, 22, 23, 59, 23, 60, 59, 0, 24, 23, 23, 24, 60, 24, 61, 60, 0, 25, 24, 24, 25, 61, 25, 62, 61, 0, 26, 25, 25, 26, 62, 26, 63, 62, 0, 27, 26, 26, 27, 63, 27, 64, 63, 0, 28, 27, 27, 28, 64, 28, 65, 64, 0, 29, 28, 28, 29, 65, 29, 66, 65, 0, 30, 29, 29, 30, 66, 30, 67, 66, 0, 31, 30, 30, 31, 67, 31, 68, 67, 0, 32, 31, 31, 32, 68, 32, 69, 68, 0, 33, 32, 32, 33, 69, 33, 70, 69, 0, 34, 33, 33, 34, 70, 34, 71, 70, 0, 35, 34, 34, 35, 71, 35, 72, 71, 0, 36, 35, 35, 36, 72, 36, 73, 72, 0, 1, 36, 36, 1, 73, 1, 38, 73, 37, 38, 39, 37, 39, 40, 37, 40, 41, 37, 41, 42, 37, 42, 43, 37, 43, 44, 37, 44, 45, 37, 45, 46, 37, 46, 47, 37, 47, 48, 37, 48, 49, 37, 49, 50, 37, 50, 51, 37, 51, 52, 37, 52, 53, 37, 53, 54, 37, 54, 55, 37, 55, 56, 37, 56, 57, 37, 57, 58, 37, 58, 59, 37, 59, 60, 37, 60, 61, 37, 61, 62, 37, 62, 63, 37, 63, 64, 37, 64, 65, 37, 65, 66, 37, 66, 67, 37, 67, 68, 37, 68, 69, 37, 69, 70, 37, 70, 71, 37, 71, 72, 37, 72, 73, 37, 73, 38], [ 0, 2, 1, 1, 2, 38, 2, 39, 38, 0, 3, 2, 2, 3, 39, 3, 40, 39, 0, 4, 3, 3, 4, 40, 4, 41, 40, 0, 5, 4, 4, 5, 41, 5, 42, 41, 0, 6, 5, 5, 6, 42, 6, 43, 42, 0, 7, 6, 6, 7, 43, 7, 44, 43, 0, 8, 7, 7, 8, 44, 8, 45, 44, 0, 9, 8, 8, 9, 45, 9, 46, 45, 0, 10, 9, 9, 10, 46, 10, 47, 46, 0, 11, 10, 10, 11, 47, 11, 48, 47, 0, 12, 11, 11, 12, 48, 12, 49, 48, 0, 13, 12, 12, 13, 49, 13, 50, 49, 0, 14, 13, 13, 14, 50, 14, 51, 50, 0, 15, 14, 14, 15, 51, 15, 52, 51, 0, 16, 15, 15, 16, 52, 16, 53, 52, 0, 17, 16, 16, 17, 53, 17, 54, 53, 0, 18, 17, 17, 18, 54, 18, 55, 54, 0, 19, 18, 18, 19, 55, 19, 56, 55, 0, 20, 19, 19, 20, 56, 20, 57, 56, 0, 21, 20, 20, 21, 57, 21, 58, 57, 0, 22, 21, 21, 22, 58, 22, 59, 58, 0, 23, 22, 22, 23, 59, 23, 60, 59, 0, 24, 23, 23, 24, 60, 24, 61, 60, 0, 25, 24, 24, 25, 61, 25, 62, 61, 0, 26, 25, 25, 26, 62, 26, 63, 62, 0, 27, 26, 26, 27, 63, 27, 64, 63, 0, 28, 27, 27, 28, 64, 28, 65, 64, 0, 29, 28, 28, 29, 65, 29, 66, 65, 0, 30, 29, 29, 30, 66, 30, 67, 66, 0, 31, 30, 30, 31, 67, 31, 68, 67, 0, 32, 31, 31, 32, 68, 32, 69, 68, 0, 33, 32, 32, 33, 69, 33, 70, 69, 0, 34, 33, 33, 34, 70, 34, 71, 70, 0, 35, 34, 34, 35, 71, 35, 72, 71, 0, 36, 35, 35, 36, 72, 36, 73, 72, 0, 1, 36, 36, 1, 73, 1, 38, 73, 37, 38, 39, 37, 39, 40, 37, 40, 41, 37, 41, 42, 37, 42, 43, 37, 43, 44, 37, 44, 45, 37, 45, 46, 37, 46, 47, 37, 47, 48, 37, 48, 49, 37, 49, 50, 37, 50, 51, 37, 51, 52, 37, 52, 53, 37, 53, 54, 37, 54, 55, 37, 55, 56, 37, 56, 57, 37, 57, 58, 37, 58, 59, 37, 59, 60, 37, 60, 61, 37, 61, 62, 37, 62, 63, 37, 63, 64, 37, 64, 65, 37, 65, 66, 37, 66, 67, 37, 67, 68, 37, 68, 69, 37, 69, 70, 37, 70, 71, 37, 71, 72, 37, 72, 73, 37, 73, 38], [ 2, 1, 0, 2, 4, 3, 2, 5, 4, 2, 0, 49, 2, 49, 5, 49, 48, 5, 5, 48, 6, 27, 26, 25, 27, 29, 28, 27, 30, 29, 27, 25, 24, 27, 24, 30, 24, 23, 30, 30, 23, 31, 52, 50, 51, 52, 53, 54, 52, 54, 55, 52, 99, 50, 52, 55, 99, 99, 55, 98, 55, 56, 98, 77, 75, 76, 77, 78, 79, 77, 79, 80, 77, 74, 75, 77, 80, 74, 74, 80, 73, 80, 81, 73, 6, 48, 7, 48, 47, 7, 56, 57, 98, 98, 57, 97, 7, 47, 8, 47, 46, 8, 57, 58, 97, 97, 58, 96, 8, 46, 9, 46, 45, 9, 58, 59, 96, 96, 59, 95, 9, 45, 10, 45, 44, 10, 59, 60, 95, 95, 60, 94, 10, 44, 11, 44, 43, 11, 60, 61, 94, 94, 61, 93, 11, 43, 12, 43, 42, 12, 61, 62, 93, 93, 62, 92, 12, 42, 13, 42, 41, 13, 62, 63, 92, 92, 63, 91, 13, 41, 14, 41, 40, 14, 63, 64, 91, 91, 64, 90, 14, 40, 15, 40, 39, 15, 64, 65, 90, 90, 65, 89, 15, 39, 16, 39, 38, 16, 65, 66, 89, 89, 66, 88, 16, 38, 17, 38, 37, 17, 66, 67, 88, 88, 67, 87, 17, 37, 18, 37, 36, 18, 67, 68, 87, 87, 68, 86, 18, 36, 19, 36, 35, 19, 68, 69, 86, 86, 69, 85, 19, 35, 20, 35, 34, 20, 69, 70, 85, 85, 70, 84, 20, 34, 21, 34, 33, 21, 70, 71, 84, 84, 71, 83, 21, 33, 22, 33, 32, 22, 71, 72, 83, 83, 72, 82, 22, 32, 23, 32, 31, 23, 72, 73, 82, 82, 73, 81, 0, 1, 51, 0, 51, 50, 1, 2, 52, 1, 52, 51, 2, 3, 53, 2, 53, 52, 3, 4, 54, 3, 54, 53, 6, 7, 57, 6, 57, 56, 7, 8, 58, 7, 58, 57, 8, 9, 59, 8, 59, 58, 9, 10, 60, 9, 60, 59, 10, 11, 61, 10, 61, 60, 11, 12, 62, 11, 62, 61, 12, 13, 63, 12, 63, 62, 13, 14, 64, 13, 64, 63, 14, 15, 65, 14, 65, 64, 15, 16, 66, 15, 66, 65, 16, 17, 67, 16, 67, 66, 17, 18, 68, 17, 68, 67, 18, 19, 69, 18, 69, 68, 19, 20, 70, 19, 70, 69, 20, 21, 71, 20, 71, 70, 21, 22, 72, 21, 72, 71, 22, 23, 73, 22, 73, 72, 25, 26, 76, 25, 76, 75, 26, 27, 77, 26, 77, 76, 27, 28, 78, 27, 78, 77, 28, 29, 79, 28, 79, 78, 31, 32, 82, 31, 82, 81, 32, 33, 83, 32, 83, 82, 33, 34, 84, 33, 84, 83, 34, 35, 85, 34, 85, 84, 35, 36, 86, 35, 86, 85, 36, 37, 87, 36, 87, 86, 37, 38, 88, 37, 88, 87, 38, 39, 89, 38, 89, 88, 39, 40, 90, 39, 90, 89, 40, 41, 91, 40, 91, 90, 41, 42, 92, 41, 92, 91, 42, 43, 93, 42, 93, 92, 43, 44, 94, 43, 94, 93, 44, 45, 95, 44, 95, 94, 45, 46, 96, 45, 96, 95, 46, 47, 97, 46, 97, 96, 47, 48, 98, 47, 98, 97, 49, 0, 48, 100, 48, 0, 100, 98, 48, 100, 50, 98, 99, 98, 50, 4, 5, 6, 4, 6, 101, 101, 6, 56, 54, 101, 56, 54, 56, 55, 24, 25, 23, 23, 25, 102, 23, 102, 73, 73, 102, 75, 74, 73, 75, 31, 29, 30, 31, 103, 29, 103, 31, 81, 79, 103, 81, 79, 81, 80], [ 0, 2, 1, 1, 2, 38, 2, 39, 38, 0, 3, 2, 2, 3, 39, 3, 40, 39, 0, 4, 3, 3, 4, 40, 4, 41, 40, 0, 5, 4, 4, 5, 41, 5, 42, 41, 0, 6, 5, 5, 6, 42, 6, 43, 42, 0, 7, 6, 6, 7, 43, 7, 44, 43, 0, 8, 7, 7, 8, 44, 8, 45, 44, 0, 9, 8, 8, 9, 45, 9, 46, 45, 0, 10, 9, 9, 10, 46, 10, 47, 46, 0, 11, 10, 10, 11, 47, 11, 48, 47, 0, 12, 11, 11, 12, 48, 12, 49, 48, 0, 13, 12, 12, 13, 49, 13, 50, 49, 0, 14, 13, 13, 14, 50, 14, 51, 50, 0, 15, 14, 14, 15, 51, 15, 52, 51, 0, 16, 15, 15, 16, 52, 16, 53, 52, 0, 17, 16, 16, 17, 53, 17, 54, 53, 0, 18, 17, 17, 18, 54, 18, 55, 54, 0, 19, 18, 18, 19, 55, 19, 56, 55, 0, 20, 19, 19, 20, 56, 20, 57, 56, 0, 21, 20, 20, 21, 57, 21, 58, 57, 0, 22, 21, 21, 22, 58, 22, 59, 58, 0, 23, 22, 22, 23, 59, 23, 60, 59, 0, 24, 23, 23, 24, 60, 24, 61, 60, 0, 25, 24, 24, 25, 61, 25, 62, 61, 0, 26, 25, 25, 26, 62, 26, 63, 62, 0, 27, 26, 26, 27, 63, 27, 64, 63, 0, 28, 27, 27, 28, 64, 28, 65, 64, 0, 29, 28, 28, 29, 65, 29, 66, 65, 0, 30, 29, 29, 30, 66, 30, 67, 66, 0, 31, 30, 30, 31, 67, 31, 68, 67, 0, 32, 31, 31, 32, 68, 32, 69, 68, 0, 33, 32, 32, 33, 69, 33, 70, 69, 0, 34, 33, 33, 34, 70, 34, 71, 70, 0, 35, 34, 34, 35, 71, 35, 72, 71, 0, 36, 35, 35, 36, 72, 36, 73, 72, 0, 1, 36, 36, 1, 73, 1, 38, 73, 37, 38, 39, 37, 39, 40, 37, 40, 41, 37, 41, 42, 37, 42, 43, 37, 43, 44, 37, 44, 45, 37, 45, 46, 37, 46, 47, 37, 47, 48, 37, 48, 49, 37, 49, 50, 37, 50, 51, 37, 51, 52, 37, 52, 53, 37, 53, 54, 37, 54, 55, 37, 55, 56, 37, 56, 57, 37, 57, 58, 37, 58, 59, 37, 59, 60, 37, 60, 61, 37, 61, 62, 37, 62, 63, 37, 63, 64, 37, 64, 65, 37, 65, 66, 37, 66, 67, 37, 67, 68, 37, 68, 69, 37, 69, 70, 37, 70, 71, 37, 71, 72, 37, 72, 73, 37, 73, 38], [ 0, 2, 1, 1, 2, 38, 2, 39, 38, 0, 3, 2, 2, 3, 39, 3, 40, 39, 0, 4, 3, 3, 4, 40, 4, 41, 40, 0, 5, 4, 4, 5, 41, 5, 42, 41, 0, 6, 5, 5, 6, 42, 6, 43, 42, 0, 7, 6, 6, 7, 43, 7, 44, 43, 0, 8, 7, 7, 8, 44, 8, 45, 44, 0, 9, 8, 8, 9, 45, 9, 46, 45, 0, 10, 9, 9, 10, 46, 10, 47, 46, 0, 11, 10, 10, 11, 47, 11, 48, 47, 0, 12, 11, 11, 12, 48, 12, 49, 48, 0, 13, 12, 12, 13, 49, 13, 50, 49, 0, 14, 13, 13, 14, 50, 14, 51, 50, 0, 15, 14, 14, 15, 51, 15, 52, 51, 0, 16, 15, 15, 16, 52, 16, 53, 52, 0, 17, 16, 16, 17, 53, 17, 54, 53, 0, 18, 17, 17, 18, 54, 18, 55, 54, 0, 19, 18, 18, 19, 55, 19, 56, 55, 0, 20, 19, 19, 20, 56, 20, 57, 56, 0, 21, 20, 20, 21, 57, 21, 58, 57, 0, 22, 21, 21, 22, 58, 22, 59, 58, 0, 23, 22, 22, 23, 59, 23, 60, 59, 0, 24, 23, 23, 24, 60, 24, 61, 60, 0, 25, 24, 24, 25, 61, 25, 62, 61, 0, 26, 25, 25, 26, 62, 26, 63, 62, 0, 27, 26, 26, 27, 63, 27, 64, 63, 0, 28, 27, 27, 28, 64, 28, 65, 64, 0, 29, 28, 28, 29, 65, 29, 66, 65, 0, 30, 29, 29, 30, 66, 30, 67, 66, 0, 31, 30, 30, 31, 67, 31, 68, 67, 0, 32, 31, 31, 32, 68, 32, 69, 68, 0, 33, 32, 32, 33, 69, 33, 70, 69, 0, 34, 33, 33, 34, 70, 34, 71, 70, 0, 35, 34, 34, 35, 71, 35, 72, 71, 0, 36, 35, 35, 36, 72, 36, 73, 72, 0, 1, 36, 36, 1, 73, 1, 38, 73, 37, 38, 39, 37, 39, 40, 37, 40, 41, 37, 41, 42, 37, 42, 43, 37, 43, 44, 37, 44, 45, 37, 45, 46, 37, 46, 47, 37, 47, 48, 37, 48, 49, 37, 49, 50, 37, 50, 51, 37, 51, 52, 37, 52, 53, 37, 53, 54, 37, 54, 55, 37, 55, 56, 37, 56, 57, 37, 57, 58, 37, 58, 59, 37, 59, 60, 37, 60, 61, 37, 61, 62, 37, 62, 63, 37, 63, 64, 37, 64, 65, 37, 65, 66, 37, 66, 67, 37, 67, 68, 37, 68, 69, 37, 69, 70, 37, 70, 71, 37, 71, 72, 37, 72, 73, 37, 73, 38], [ 0, 2, 1, 1, 2, 38, 2, 39, 38, 0, 3, 2, 2, 3, 39, 3, 40, 39, 0, 4, 3, 3, 4, 40, 4, 41, 40, 0, 5, 4, 4, 5, 41, 5, 42, 41, 0, 6, 5, 5, 6, 42, 6, 43, 42, 0, 7, 6, 6, 7, 43, 7, 44, 43, 0, 8, 7, 7, 8, 44, 8, 45, 44, 0, 9, 8, 8, 9, 45, 9, 46, 45, 0, 10, 9, 9, 10, 46, 10, 47, 46, 0, 11, 10, 10, 11, 47, 11, 48, 47, 0, 12, 11, 11, 12, 48, 12, 49, 48, 0, 13, 12, 12, 13, 49, 13, 50, 49, 0, 14, 13, 13, 14, 50, 14, 51, 50, 0, 15, 14, 14, 15, 51, 15, 52, 51, 0, 16, 15, 15, 16, 52, 16, 53, 52, 0, 17, 16, 16, 17, 53, 17, 54, 53, 0, 18, 17, 17, 18, 54, 18, 55, 54, 0, 19, 18, 18, 19, 55, 19, 56, 55, 0, 20, 19, 19, 20, 56, 20, 57, 56, 0, 21, 20, 20, 21, 57, 21, 58, 57, 0, 22, 21, 21, 22, 58, 22, 59, 58, 0, 23, 22, 22, 23, 59, 23, 60, 59, 0, 24, 23, 23, 24, 60, 24, 61, 60, 0, 25, 24, 24, 25, 61, 25, 62, 61, 0, 26, 25, 25, 26, 62, 26, 63, 62, 0, 27, 26, 26, 27, 63, 27, 64, 63, 0, 28, 27, 27, 28, 64, 28, 65, 64, 0, 29, 28, 28, 29, 65, 29, 66, 65, 0, 30, 29, 29, 30, 66, 30, 67, 66, 0, 31, 30, 30, 31, 67, 31, 68, 67, 0, 32, 31, 31, 32, 68, 32, 69, 68, 0, 33, 32, 32, 33, 69, 33, 70, 69, 0, 34, 33, 33, 34, 70, 34, 71, 70, 0, 35, 34, 34, 35, 71, 35, 72, 71, 0, 36, 35, 35, 36, 72, 36, 73, 72, 0, 1, 36, 36, 1, 73, 1, 38, 73, 37, 38, 39, 37, 39, 40, 37, 40, 41, 37, 41, 42, 37, 42, 43, 37, 43, 44, 37, 44, 45, 37, 45, 46, 37, 46, 47, 37, 47, 48, 37, 48, 49, 37, 49, 50, 37, 50, 51, 37, 51, 52, 37, 52, 53, 37, 53, 54, 37, 54, 55, 37, 55, 56, 37, 56, 57, 37, 57, 58, 37, 58, 59, 37, 59, 60, 37, 60, 61, 37, 61, 62, 37, 62, 63, 37, 63, 64, 37, 64, 65, 37, 65, 66, 37, 66, 67, 37, 67, 68, 37, 68, 69, 37, 69, 70, 37, 70, 71, 37, 71, 72, 37, 72, 73, 37, 73, 38]]
decoration = [[ 0.000000, 0.000000, 32.632000, 0.000000, 0.000000, 33.000000], [ 0.000000, 0.000000, 33.202000, 0.000000, 0.000000, 33.570000], [ 0.000000, 0.000000, 33.772000, 0.000000, 0.000000, 34.140000], [ 0.000000, 0.000000, 34.452500, 0.000000, 0.000000, 34.462500], [ 0.000000, 0.000000, 34.452500, 0.000000, 0.000000, 34.462500], [ 0.000000, 0.000000, 36.192100, 0.000000, 0.000000, 36.202100], [ 0.000000, 0.000000, 36.192100, 0.000000, 0.000000, 36.202100], [ 0.000000, 0.000000, 36.287600, 0.000000, 0.000000, 36.297600], [ 0.000000, 0.000000, 36.287600, 0.000000, 0.000000, 36.297600], [ 0.000000, 0.000000, 37.366700, 0.000000, 0.000000, 37.376700], [ 0.000000, 0.000000, 37.366700, 0.000000, 0.000000, 37.376700], [ 0.000000, 0.000000, 37.650200, 0.000000, 0.000000, 37.970200], [ 0.000000, 0.000000, 38.280200, 0.000000, 0.000000, 38.600200], [ -0.393380, 0.000000, 39.036200, 0.393380, 0.000000, 39.036200, -0.924534, 0.000000, 40.021871, -0.491878, 0.000000, 40.678985, 0.000000, 0.000000, 39.036200, 0.000000, 0.000000, 39.036200, -0.708206, 0.000000, 40.350428, -0.708206, 0.000000, 40.350428], [ -1.036641, 0.000000, 40.566676, -1.044993, 0.000000, 40.572175], [ -1.036641, 0.000000, 40.566676, -1.044993, 0.000000, 40.572175], [ -1.263569, 0.000000, 40.716089, -1.530839, 0.000000, 40.892064], [ -1.814812, 0.000000, 41.079038, -2.082081, 0.000000, 41.255013], [ -2.220642, 0.000000, 41.817405, -2.653454, -0.000000, 41.160054, -2.751988, 0.000000, 42.802973, -3.539031, 0.000000, 42.802667, -2.437048, -0.000000, 41.488730, -2.437048, -0.000000, 41.488730, -3.145509, 0.000000, 42.802820, -3.145509, 0.000000, 42.802820], [ -3.145662, 0.000000, 43.196042, -3.145666, 0.000000, 43.206042], [ -3.145662, 0.000000, 43.196042, -3.145666, 0.000000, 43.206042], [ -3.145847, 0.000000, 43.672542, -3.145972, 0.000000, 43.992542], [ -3.146026, 0.000000, 44.132542, -3.146151, 0.000000, 44.452542], [ -3.146205, 0.000000, 44.592542, -3.146329, 0.000000, 44.912542], [ -2.757046, 0.000000, 46.443994, -3.536802, 0.000000, 46.443691, -1.617390, 0.000000, 47.613697, -1.637375, 0.000000, 48.393197, -3.146924, 0.000000, 46.443842, -3.146924, 0.000000, 46.443842, -1.627382, 0.000000, 48.003447, -1.627382, 0.000000, 48.003447], [ -1.250718, 0.000000, 48.013104, -1.240722, 0.000000, 48.013360], [ -1.250718, 0.000000, 48.013104, -1.240722, 0.000000, 48.013360], [ 0.264581, 0.000000, 48.051953, 0.274578, 0.000000, 48.052209]]

color = [2, 2, 2, 0, 0, 0, 0, 0, 0, 0, 0, 2, 2, 1, 0, 0, 2, 2, 1, 0, 0, 2, 2, 2, 1, 0, 0, 0]

index_base64 = 'eJylVlFv2zYQfs+vYGUMUACLcjIUKBLbgOMGWIBkCdJgQDfsgSZPFguJ9CjKrjf0v/dIKrYk28G2CrAl6u6+++6Od9T43cfH+cvnp1uS27KYno3DjeA1zoGJ8OiXJViGWnaVwF+1XE+iuVYWlE1etiuICA+rSWThq00dzDXhOTMV2Elts+RDRNLp2R7OSlvA9IYttoVW9EtFKlauCkAcAeM0SPfa75KEtHSTpCWruJErSyrDJ5Hjd5WmXCi6COpfKsp1meZMCVpKJTMJIrmglwgTTcdpMP6vaJwp5Zn8fwgtS/1DAItdOo5j2G07g+5yRRmShRZb8k9H4C69BpMVenNFcikEqOsDjY0UNr8iF6PRT4fCHOQyt6ekJTNLqa7I6FC0YkJItTyQfTvrLAcGlAAzZ2rNqiP0f4Cc1TXPE8at1EgR6wp9Jq9PmOZ9VsdpaJCxS2jzigd6UkyiNl9XoCBq9PrVWjNDGtsJEZrXJbYSXYK9LcA93mzvRNyFPL/uWIPCBANaK9iQm9nN5/vHX+mtfxkH5CGxpgY063o1wCx84uCNs1r5NJD4vJdjp1s1Wm0X3jIO3tvY7krTma8tYaRwBTjIvAP1kh7ok5bK3jtBHD2WSkbDjvg34Fabn+PRkCRY0iEZnQ8DufPD6qap90CFzLK66tOf68Ij0Q+I0vwdAQkQ1Qp4XSDnExgXzjz8vUFkaXSthDc6AYQILqaT2VRkZjh51hYrR+asBMOOppZ7Uc8JmgbLYBhH4R452uj2gdmcPt0NyXtc9BJOfwej4xO5Ds4os5bx3B0NRhe7nZexojpWnMZokwMUTwa4rNzmm5DLUS94F065LaHKUToYfJy9zAaD615+yEsOJJOmsjhVDEJb8G1FFkCw9AI7nbgaymxLNrnkOfGA+FaWK20sJb+AAbKBZk1YUXgVqDqeOnv/XjNsSnrnLR5QOY4wl+4X+DbZGrZ6C8vx4FH7TdaE8QkssRiKxamJjzrzq6aayNbuwgw8MTLn6VR2Gxi/DYLfP0Z/9vJ7kGvcIEay4qLf7hbPUWbEQyOP/XFfG7iITrfgPrQdLN1IA5mrETpwc+nQqs2Wvhqi9g7jX5pQKGVVyTW81XL00u3/8NeL4Fu/EQ1gwCpEu1fF86ozWF+HZWvAxm2kMDKpqdWzH+v3Wq/i0/PXw9FwAsQthh12uHeeoZJ/w+7NRiqhNxTP2Ns1uJla4UcaIuBZ4vSi4Rsj/5Wh1+z7dPf9V8c4DYfg2H/4Tb8DrvXVsg=='

def decodeVertices():
    vertices_binary = zlib.decompress(base64.b64decode(vertices_base64))
    k = 0
    for i in range(len(numVertices)):
        current = []
        for j in range(numVertices[i] * 3):
            current.append(float(struct.unpack('=d', vertices_binary[k:k+8])[0]))
            k += 8
        vertices.append(current)

def dot(a, b):
    return sum(a[i]*b[i] for i in range(len(a)))

def normalize(a):
    length = math.sqrt(dot(a, a))
    for i in range(3):
        a[i]/=length
    return a

def distance(a, b):
    c = [0] * 2
    for i in range(len(a)):
        c[i] = a[i] - b[i]
    return math.sqrt(dot(c, c))

def cross(a, b):
    c = [0, 0, 0]
    c[0] = a[1]*b[2] - a[2]*b[1]
    c[1] = a[2]*b[0] - a[0]*b[2]
    c[2] = a[0]*b[1] - a[1]*b[0]
    return c

class Quaternion:
    def __init__(self, values):
        self.R = values

    def __str__(self):
        return str(self.R)

    def __mul__(self, other):
        imagThis = self.imag()
        imagOther = other.imag()
        cro = cross(imagThis, imagOther)

        ret = ([self.R[0] * other.R[0] - dot(imagThis, imagOther)] +
               [self.R[0] * imagOther[0] + other.R[0] * imagThis[0] + cro[0],
                self.R[0] * imagOther[1] + other.R[0] * imagThis[1] + cro[1],
                self.R[0] * imagOther[2] + other.R[0] * imagThis[2] + cro[2]])

        return Quaternion(ret)

    def imag(self):
        return self.R[1:]

    def conjugate(self):
        ret = [0] * 4
        ret[0] = self.R[0]
        ret[1] = -self.R[1]
        ret[2] = -self.R[2]
        ret[3] = -self.R[3]

        return Quaternion(ret)

    def inverse(self):
        return self.conjugate()

    def rotate(self, vec3D):
        quat = Quaternion([0.0] + vec3D)
        conj = self.conjugate()

        ret = self * (quat * conj)
        return ret.R[1:]

def getQuaternion(u, ref):
    u = normalize(u)
    ref = normalize(ref)
    axis = cross(u, ref)
    normAxis = math.sqrt(dot(axis, axis))

    if normAxis < 1e-12:
        if math.fabs(dot(u, ref) - 1.0) < 1e-12:
            return Quaternion([1.0, 0.0, 0.0, 0.0])

        return Quaternion([0.0, 1.0, 0.0, 0.0])

    axis = normalize(axis)
    cosAngle = math.sqrt(0.5 * (1 + dot(u, ref)))
    sinAngle = math.sqrt(1 - cosAngle * cosAngle)

    return Quaternion([cosAngle, sinAngle * axis[0], sinAngle * axis[1], sinAngle * axis[2]])

def exportVTK():
    vertices_str = ""
    triangles_str = ""
    color_str = ""
    cellTypes_str = ""
    startIdx = 0
    vertexCounter = 0
    cellCounter = 0
    lookup_table = []
    lookup_table.append([0.5, 0.5, 0.5, 0.5])
    lookup_table.append([1.0, 0.847, 0.0, 1.0])
    lookup_table.append([1.0, 0.0, 0.0, 1.0])
    lookup_table.append([0.537, 0.745, 0.525, 1.0])
    lookup_table.append([0.5, 0.5, 0.0, 1.0])
    lookup_table.append([1.0, 0.541, 0.0, 1.0])
    lookup_table.append([0.0, 0.0, 1.0, 1.0])

    decodeVertices()

    for i in range(len(vertices)):
        for j in range(0, len(vertices[i]), 3):
            vertices_str += ("%f %f %f\n" %(vertices[i][j], vertices[i][j+1], vertices[i][j+2]))
            vertexCounter += 1

        for j in range(0, len(triangles[i]), 3):
            triangles_str += ("3 %d %d %d\n" % (triangles[i][j] + startIdx, triangles[i][j+1] + startIdx, triangles[i][j+2] + startIdx))
            cellTypes_str += "5\n"
            tmp_color = lookup_table[color[i]]
            color_str += ("%f %f %f %f\n" % (tmp_color[0], tmp_color[1], tmp_color[2], tmp_color[3]))
            cellCounter += 1
        startIdx = vertexCounter

    fh = open('70MeV_Gantry2_ElementPositions.vtk','w')
    fh.write("# vtk DataFile Version 2.0\n")
    fh.write("test\nASCII\n\n")
    fh.write("DATASET UNSTRUCTURED_GRID\n")
    fh.write("POINTS " + str(vertexCounter) + " float\n")
    fh.write(vertices_str)
    fh.write("CELLS " + str(cellCounter) + " " + str(cellCounter * 4) + "\n")
    fh.write(triangles_str)
    fh.write("CELL_TYPES " + str(cellCounter) + "\n")
    fh.write(cellTypes_str)
    fh.write("CELL_DATA " + str(cellCounter) + "\n")
    fh.write("COLOR_SCALARS type 4\n")
    fh.write(color_str + "\n")
    fh.close()

def getNormal(tri_vertices):
    vec1 = [tri_vertices[1][0] - tri_vertices[0][0],
            tri_vertices[1][1] - tri_vertices[0][1],
            tri_vertices[1][2] - tri_vertices[0][2]]
    vec2 = [tri_vertices[2][0] - tri_vertices[0][0],
            tri_vertices[2][1] - tri_vertices[0][1],
            tri_vertices[2][2] - tri_vertices[0][2]]
    return normalize(cross(vec1,vec2))

def exportWeb(bgcolor):
    lookup_table = []
    lookup_table.append([0.5, 0.5, 0.5])
    lookup_table.append([1.0, 0.847, 0.0])
    lookup_table.append([1.0, 0.0, 0.0])
    lookup_table.append([0.537, 0.745, 0.525])
    lookup_table.append([0.5, 0.5, 0.0])
    lookup_table.append([1.0, 0.541, 0.0])
    lookup_table.append([0.0, 0.0, 1.0])

    decodeVertices()

    mesh = "'data:"
    mesh += "{"
    mesh += "\"autoClear\":true,"
    mesh += "\"clearColor\":[0.0000,0.0000,0.0000],"
    mesh += "\"ambientColor\":[0.0000,0.0000,0.0000],"
    mesh += "\"gravity\":[0.0000,-9.8100,0.0000],"
    mesh += "\"cameras\":["
    mesh += "{"
    mesh += "\"name\":\"Camera\","
    mesh += "\"id\":\"Camera\","
    mesh += "\"position\":[21.7936,2.2312,-85.7292],"
    mesh += "\"rotation\":[0.0432,-0.1766,-0.0668],"
    mesh += "\"fov\":0.8578,"
    mesh += "\"minZ\":10.0000,"
    mesh += "\"maxZ\":10000.0000,"
    mesh += "\"speed\":1.0000,"
    mesh += "\"inertia\":0.9000,"
    mesh += "\"checkCollisions\":false,"
    mesh += "\"applyGravity\":false,"
    mesh += "\"ellipsoid\":[0.2000,0.9000,0.2000]"
    mesh += "}],"
    mesh += "\"activeCamera\":\"Camera\","
    mesh += "\"lights\":["
    mesh += "{"
    mesh += "\"name\":\"Lamp\","
    mesh += "\"id\":\"Lamp\","
    mesh += "\"type\":0.0000,"
    mesh += "\"position\":[4.0762,34.9321,-63.5788],"
    mesh += "\"intensity\":1.0000,"
    mesh += "\"diffuse\":[1.0000,1.0000,1.0000],"
    mesh += "\"specular\":[1.0000,1.0000,1.0000]"
    mesh += "}],"
    mesh += "\"materials\":[],"
    mesh += "\"meshes\":["

    for i in range(len(triangles)):
        vertex_list = []
        indices_list = []
        normals_list = []
        color_list = []
        for j in range(0, len(triangles[i]), 3):
            tri_vertices = []
            idcs = triangles[i][j:j + 3]
            tri_vertices.append(vertices[i][3 * idcs[0]:3 * (idcs[0] + 1)])
            tri_vertices.append(vertices[i][3 * idcs[1]:3 * (idcs[1] + 1)])
            tri_vertices.append(vertices[i][3 * idcs[2]:3 * (idcs[2] + 1)])
            indices_list.append(','.join(str(n) for n in range(len(vertex_list),len(vertex_list) + 3)))
            # left hand order!
            vertex_list.append(','.join("%.5f" % (round(n,5) + 0.0) for n in tri_vertices[0]))
            vertex_list.append(','.join("%.5f" % (round(n,5) + 0.0) for n in tri_vertices[2]))
            vertex_list.append(','.join("%.5f" % (round(n,5) + 0.0) for n in tri_vertices[1]))
            normal = getNormal(tri_vertices)
            normals_list.append(','.join("%.5f" % (round(n,5) + 0.0) for n in normal * 3))
            color_list.append(','.join([str(n) for n in lookup_table[color[i]]] * 3))
        mesh += "{"
        mesh += "\"name\":\"element_" + str(i) + "\","
        mesh += "\"id\":\"element_" + str(i) + "\","
        mesh += "\"position\":[0.0,0.0,0.0],"
        mesh += "\"rotation\":[0.0,0.0,0.0],"
        mesh += "\"scaling\":[1.0,1.0,1.0],"
        mesh += "\"isVisible\":true,"
        mesh += "\"isEnabled\":true,"
        mesh += "\"useFlatShading\":false,"
        mesh += "\"checkCollisions\":false,"
        mesh += "\"billboardMode\":0,"
        mesh += "\"receiveShadows\":false,"
        mesh += "\"positions\":[" + ','.join(vertex_list) + "],"
        mesh += "\"normals\":[" + ','.join(normals_list) + "],"
        mesh += "\"indices\":[" + ','.join(indices_list) + "],"
        mesh += "\"colors\":[" + ','.join(color_list) + "],"
        mesh += "\"subMeshes\":["
        mesh += "{"
        mesh += "\"materialIndex\":0,"
        mesh += " \"verticesStart\":0,"
        mesh += " \"verticesCount\":" + str(len(triangles[i])) + ","
        mesh += " \"indexStart\":0,"
        mesh += " \"indexCount\":" + str(len(triangles[i])) + ""
        mesh += "}]"
        mesh += "}"
        mesh += ","

        del normals_list[:]
        del vertex_list[:]
        del color_list[:]
        del indices_list[:]

    mesh = mesh[:-1] + "]"
    mesh += "}'"
    index_compressed = base64.b64decode(index_base64)
    index = str(zlib.decompress(index_compressed))
    if (len(bgcolor) == 3):
        mesh += ";\n            "
        mesh += "scene.clearColor = new BABYLON.Color3(%f, %f, %f)" % (bgcolor[0], bgcolor[1], bgcolor[2])

    index = index.replace('##DATA##', mesh)
    fh = open('70MeV_Gantry2_ElementPositions.html','w')
    fh.write(index)
    fh.close()

def computeMinAngle(idx, curAngle, positions, connections, check):
    matrix = [[-math.cos(curAngle), -math.sin(curAngle)], [math.sin(curAngle), -math.cos(curAngle)]]

    minAngle = 2 * math.pi
    nextIdx = -1

    for j in connections[idx]:
        direction = [positions[j][0] - positions[idx][0],
                     positions[j][1] - positions[idx][1]]
        direction = [dot(matrix[0],direction), dot(matrix[1],direction)]

        if math.fabs(dot([1.0, 0.0], direction) / distance(positions[j], positions[idx]) - 1.0) < 1e-6: continue

        angle = math.fmod(math.atan2(direction[1],direction[0]) + 2 * math.pi, 2 * math.pi)

        if angle < minAngle:
            nextIdx = j
            minAngle = angle
        if angle == minAngle and check:
            dire =  [positions[j][0] - positions[idx][0],
                     positions[j][1] - positions[idx][1]]
            minA0 = math.atan2(dire[1], dire[0])
            minA1 = computeMinAngle(nextIdx, minA0, positions, connections, False)
            minA2 = computeMinAngle(j, minA0, positions, connections, False)
            if minA2[1] < minA2[1]:
                nextIdx = j

    if nextIdx == -1:
        nextIdx = connections[idx][0]

    return (nextIdx, minAngle)

def squashVertices(positionDict, connections):
    removedItems = []
    idxChanges = positionDict.keys()
    for i in positionDict.keys():
        if i in removedItems:
            continue
        for j in positionDict.keys():
            if j in removedItems or j == i:
                continue

            if distance(positionDict[i], positionDict[j]) < 1e-6:
                connections[j] = list(set(connections[j]))
                if i in connections[j]:
                    connections[j].remove(i)
                if j in connections[i]:
                    connections[i].remove(j)

                connections[i].extend(connections[j])
                connections[i] = list(set(connections[i]))
                connections[i].sort()

                for k in connections.keys():
                    if k == i: continue
                    if j in connections[k]:
                        idx = connections[k].index(j)
                        connections[k][idx] = i

                idxChanges[j] = i
                removedItems.append(j)

    for i in removedItems:
        del positionDict[i]
        del connections[i]

    for i in connections.keys():
        connections[i] = list(set(connections[i]))
        connections[i].sort()

    return idxChanges

def computeLineEquations(positions, connections):
    cons = set()
    for i in connections.keys():
        for j in connections[i]:
            cons.add((min(i, j), max(i, j)))

    lineEquations = {}
    for item in cons:
        a = (positions[item[1]][1] - positions[item[0]][1])
        b = -(positions[item[1]][0] - positions[item[0]][0])

        xm = 0.5 * (positions[item[0]][0] + positions[item[1]][0])
        ym = 0.5 * (positions[item[0]][1] + positions[item[1]][1])
        c = -(a * xm +  b * ym)

        lineEquations[item] = (a, b, c)

    return lineEquations

def checkPossibleSegmentIntersection(segment, positions, lineEquations, minAngle, lastIntersection):
    item1 = (min(segment), max(segment))

    (a1, b1, c1) = (0,0,0)
    A = [0]*2
    B = A

    if segment[0] == None:
        A = lastIntersection
        B = positions[segment[1]]

        a1 = B[1] - A[1]
        b1 = -(B[0] - A[0])
        xm = 0.5 * (A[0] + B[0])
        ym = 0.5 * (A[1] + B[1])
        c1 = -(a1 * xm + b1 * ym)

    else:
        A = positions[segment[0]]
        B = positions[segment[1]]

        (a1, b1, c1) = lineEquations[item1]

        if segment[1] < segment[0]:
            (a1, b1, c1) = (-a1, -b1, -c1)

    curAngle = math.atan2(a1, -b1)
    matrix = [[-math.cos(curAngle), -math.sin(curAngle)], [math.sin(curAngle), -math.cos(curAngle)]]

    origMinAngle = minAngle

    segment1 = [B[0] - A[0], B[1] - A[1], 0.0]

    intersectingSegments = []
    distanceAB = distance(A, B)
    for item2 in lineEquations.keys():
        item = item2
        C = positions[item[0]]
        D = positions[item[1]]

        if segment[0] == None:
            if (segment[1] == item[0] or
                segment[1] == item[1]): continue
        else:
            if (item1[0] == item[0] or
                item1[1] == item[0] or
                item1[0] == item[1] or
                item1[1] == item[1]): continue

        nodes = set([item1[0],item1[1],item[0],item[1]])
        if len(nodes) < 4: continue       # share same vertex

        (a2, b2, c2) = lineEquations[item]

        segment2 = [C[0] - A[0], C[1] - A[1], 0.0]
        segment3 = [D[0] - A[0], D[1] - A[1], 0.0]

        # check that C and D aren't on the same side of AB
        if cross(segment1, segment2)[2] * cross(segment1, segment3)[2] > 0.0: continue

        if cross(segment1, segment2)[2] < 0.0 or cross(segment1, segment3)[2] > 0.0:
            (C, D, a2, b2, c2) = (D, C, -a2, -b2, -c2)
            item = (item[1], item[0])

        denominator = a1 * b2 - b1 * a2
        if math.fabs(denominator) < 1e-9: continue

        px = (b1 * c2 - b2 * c1) / denominator
        py = (a2 * c1 - a1 * c2) / denominator

        distanceCD = distance(C, D)

        distanceAP = distance(A, [px, py])
        distanceBP = distance(B, [px, py])
        distanceCP = distance(C, [px, py])
        distanceDP = distance(D, [px, py])

        # check intersection is between AB and CD
        check1 = (distanceAP - 1e-6 < distanceAB and distanceBP - 1e-6 < distanceAB)
        check2 = (distanceCP - 1e-6 < distanceCD and distanceDP - 1e-6 < distanceCD)
        if not check1 or not check2: continue

        if math.fabs(dot(segment1, [D[0] - C[0], D[1] - C[1], 0.0]) / (distanceAB * distanceCD) + 1.0) < 1e-9: continue

        if ((distanceAP < 1e-6) or
            (distanceBP < 1e-6 and distanceCP < 1e-6) or
            (distanceDP < 1e-6)):
            continue

        direction = [D[0] - C[0], D[1] - C[1]]
        direction = [dot(matrix[0], direction), dot(matrix[1], direction)]
        angle = math.fmod(math.atan2(direction[1], direction[0]) + 2 * math.pi, 2 * math.pi)

        newSegment = ([px, py], item[1], distanceAP, angle)

        if distanceCP < 1e-6 and angle > origMinAngle: continue
        if distanceBP > 1e-6 and angle > math.pi: continue

        if len(intersectingSegments) == 0:
            intersectingSegments.append(newSegment)
            minAngle = angle
        else:
            if intersectingSegments[0][2] - 1e-9 > distanceAP:
                intersectingSegments[0] = newSegment
                minAngle = angle
            elif intersectingSegments[0][2] + 1e-9 > distanceAP and angle < minAngle:
                intersectingSegments[0] = newSegment
                minAngle = angle

    return intersectingSegments

def projectToPlane(normal):
    fh = open('70MeV_Gantry2_ElementPositions.gpl','w')
    normal = normalize(normal)
    ori = getQuaternion(normal, [0, 0, 1])

    left2Right = [0, 0, 1]
    if math.fabs(math.fabs(dot(normal, left2Right)) - 1) < 1e-9:
        left2Right = [1, 0, 0]
    rotL2R = ori.rotate(left2Right)
    deviation = math.atan2(rotL2R[1], rotL2R[0])
    rotBack = Quaternion([math.cos(0.5 * deviation), 0, 0, -math.sin(0.5 * deviation)])
    ori = rotBack * ori

    decodeVertices()

    for i in range(len(vertices)):
        positions = {}
        connections = {}
        for j in range(0, len(vertices[i]), 3):
            nextPos3D = ori.rotate(vertices[i][j:j+3])
            nextPos2D = nextPos3D[0:2]
            positions[j/3] = nextPos2D

        if len(positions) == 0:
            continue

        idx = 0
        maxX = positions[0][0]
        for j in positions.keys():
            if positions[j][0] > maxX:
                maxX = positions[j][0]
                idx = j
            if positions[j][0] == maxX and positions[j][1] > positions[idx][1]:
                idx = j

        for j in range(0, len(triangles[i]), 3):
            for k in range(0, 3):
                vertIdx = triangles[i][j + k]
                if not vertIdx in connections:
                    connections[vertIdx] = []
                for l in range(1, 3):
                    connections[vertIdx].append(triangles[i][j + ((k + l) % 3)])

        numConnections = 0
        for j in connections.keys():
            connections[j] = list(set(connections[j]))
            numConnections += len(connections[j])

        numConnections /= 2
        idChanges = squashVertices(positions, connections)

        lineEquations = computeLineEquations(positions, connections)

        idx = idChanges[idx]
        fh.write("%.6f    %.6f    %d\n" % (positions[idx][0], positions[idx][1], idx))

        curAngle = math.pi
        origin = idx
        count = 0
        passed = []
        while (count == 0 or distance(positions[origin], positions[idx]) > 1e-9) and not count > numConnections:
            nextGen = computeMinAngle(idx, curAngle, positions, connections, False)
            nextIdx = nextGen[0]
            direction = [positions[nextIdx][0] - positions[idx][0],
                         positions[nextIdx][1] - positions[idx][1]]
            curAngle = math.atan2(direction[1], direction[0])

            intersections = checkPossibleSegmentIntersection((idx, nextIdx), positions, lineEquations, nextGen[1], [])
            if len(intersections) > 0:
                while len(intersections) > 0 and not count > numConnections:
                    fh.write("%.6f    %.6f\n" %(intersections[0][0][0], intersections[0][0][1]))
                    count += 1
    
                    idx = intersections[0][1]
                    direction = [positions[idx][0] - intersections[0][0][0],
                                 positions[idx][1] - intersections[0][0][1]]
                    curAngle = math.atan2(direction[1], direction[0])
                    nextGen = computeMinAngle(idx, curAngle, positions, connections, False)
    
                    intersections = checkPossibleSegmentIntersection((None, idx), positions, lineEquations, nextGen[1], intersections[0][0])
            else:
                idx = nextIdx
    
            fh.write("%.6f    %.6f    %d\n" % (positions[idx][0], positions[idx][1], idx))
    
            if idx in passed:
                direction1 = [positions[idx][0] - positions[passed[-1]][0],
                              positions[idx][1] - positions[passed[-1]][1]]
                direction2 = [positions[origin][0] - positions[passed[-1]][0],
                              positions[origin][1] - positions[passed[-1]][1]]
                dist1 = distance(positions[idx], positions[passed[-1]])
                dist2 = distance(positions[origin], positions[passed[-1]])
                if dist1 * dist2 > 0.0 and math.fabs(math.fabs(dot(direction1, direction2) / (dist1 * dist2)) - 1.0) > 1e-9:
                    sys.stderr.write("error: projection cycling on element id: %d, vertex id: %d\n" %(i, idx))
                break

            passed.append(idx)
            count += 1
        fh.write("\n")

        if count > numConnections:
            sys.stderr.write("error: projection cycling on element id: %d\n" % i)

        for j in range(0, len(decoration[i]), 6):
            for k in range(j, j + 6, 3):
                nextPos3D = ori.rotate(decoration[i][k:k+3])
                fh.write("%.6f    %.6f\n" % (nextPos3D[0], nextPos3D[1]))
            fh.write("\n")

    fh.close()

parser = argparse.ArgumentParser()
parser.add_argument('--export-vtk', action='store_true')
parser.add_argument('--export-web', action='store_true')
parser.add_argument('--background', nargs=3, type=float)
parser.add_argument('--project-to-plane', action='store_true')
parser.add_argument('--normal', nargs=3, type=float)
args = parser.parse_args()

if (args.export_vtk):
    exportVTK()
    sys.exit()

if (args.export_web):
    bgcolor = []
    if (args.background):
        validBackground = True
        for comp in bgcolor:
            if comp < 0.0 or comp > 1.0:
                validBackground = False
                break
        if (validBackground):
            bgcolor = args.background
    exportWeb(bgcolor)
    sys.exit()

if (args.project_to_plane):
    normal = [0.0, 1.0, 0.0]
    if (args.normal):
        normal = args.normal
    projectToPlane(normal)
    sys.exit()

parser.print_help()
