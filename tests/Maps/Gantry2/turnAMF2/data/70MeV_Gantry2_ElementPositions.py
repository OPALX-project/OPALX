import os, sys, argparse, math, base64, zlib, struct

if sys.version_info < (3,0):
    range = xrange

vertices = []
numVertices = [74, 74, 74, 74, 74, 74, 74, 74, 74, 74, 74, 74, 74, 80, 74, 74, 74, 74, 80, 74, 74, 74, 74, 74, 104, 74, 74, 74]
vertices_base64 = "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"
triangles = [[ 0, 2, 1, 1, 2, 38, 2, 39, 38, 0, 3, 2, 2, 3, 39, 3, 40, 39, 0, 4, 3, 3, 4, 40, 4, 41, 40, 0, 5, 4, 4, 5, 41, 5, 42, 41, 0, 6, 5, 5, 6, 42, 6, 43, 42, 0, 7, 6, 6, 7, 43, 7, 44, 43, 0, 8, 7, 7, 8, 44, 8, 45, 44, 0, 9, 8, 8, 9, 45, 9, 46, 45, 0, 10, 9, 9, 10, 46, 10, 47, 46, 0, 11, 10, 10, 11, 47, 11, 48, 47, 0, 12, 11, 11, 12, 48, 12, 49, 48, 0, 13, 12, 12, 13, 49, 13, 50, 49, 0, 14, 13, 13, 14, 50, 14, 51, 50, 0, 15, 14, 14, 15, 51, 15, 52, 51, 0, 16, 15, 15, 16, 52, 16, 53, 52, 0, 17, 16, 16, 17, 53, 17, 54, 53, 0, 18, 17, 17, 18, 54, 18, 55, 54, 0, 19, 18, 18, 19, 55, 19, 56, 55, 0, 20, 19, 19, 20, 56, 20, 57, 56, 0, 21, 20, 20, 21, 57, 21, 58, 57, 0, 22, 21, 21, 22, 58, 22, 59, 58, 0, 23, 22, 22, 23, 59, 23, 60, 59, 0, 24, 23, 23, 24, 60, 24, 61, 60, 0, 25, 24, 24, 25, 61, 25, 62, 61, 0, 26, 25, 25, 26, 62, 26, 63, 62, 0, 27, 26, 26, 27, 63, 27, 64, 63, 0, 28, 27, 27, 28, 64, 28, 65, 64, 0, 29, 28, 28, 29, 65, 29, 66, 65, 0, 30, 29, 29, 30, 66, 30, 67, 66, 0, 31, 30, 30, 31, 67, 31, 68, 67, 0, 32, 31, 31, 32, 68, 32, 69, 68, 0, 33, 32, 32, 33, 69, 33, 70, 69, 0, 34, 33, 33, 34, 70, 34, 71, 70, 0, 35, 34, 34, 35, 71, 35, 72, 71, 0, 36, 35, 35, 36, 72, 36, 73, 72, 0, 1, 36, 36, 1, 73, 1, 38, 73, 37, 38, 39, 37, 39, 40, 37, 40, 41, 37, 41, 42, 37, 42, 43, 37, 43, 44, 37, 44, 45, 37, 45, 46, 37, 46, 47, 37, 47, 48, 37, 48, 49, 37, 49, 50, 37, 50, 51, 37, 51, 52, 37, 52, 53, 37, 53, 54, 37, 54, 55, 37, 55, 56, 37, 56, 57, 37, 57, 58, 37, 58, 59, 37, 59, 60, 37, 60, 61, 37, 61, 62, 37, 62, 63, 37, 63, 64, 37, 64, 65, 37, 65, 66, 37, 66, 67, 37, 67, 68, 37, 68, 69, 37, 69, 70, 37, 70, 71, 37, 71, 72, 37, 72, 73, 37, 73, 38], [ 0, 2, 1, 1, 2, 38, 2, 39, 38, 0, 3, 2, 2, 3, 39, 3, 40, 39, 0, 4, 3, 3, 4, 40, 4, 41, 40, 0, 5, 4, 4, 5, 41, 5, 42, 41, 0, 6, 5, 5, 6, 42, 6, 43, 42, 0, 7, 6, 6, 7, 43, 7, 44, 43, 0, 8, 7, 7, 8, 44, 8, 45, 44, 0, 9, 8, 8, 9, 45, 9, 46, 45, 0, 10, 9, 9, 10, 46, 10, 47, 46, 0, 11, 10, 10, 11, 47, 11, 48, 47, 0, 12, 11, 11, 12, 48, 12, 49, 48, 0, 13, 12, 12, 13, 49, 13, 50, 49, 0, 14, 13, 13, 14, 50, 14, 51, 50, 0, 15, 14, 14, 15, 51, 15, 52, 51, 0, 16, 15, 15, 16, 52, 16, 53, 52, 0, 17, 16, 16, 17, 53, 17, 54, 53, 0, 18, 17, 17, 18, 54, 18, 55, 54, 0, 19, 18, 18, 19, 55, 19, 56, 55, 0, 20, 19, 19, 20, 56, 20, 57, 56, 0, 21, 20, 20, 21, 57, 21, 58, 57, 0, 22, 21, 21, 22, 58, 22, 59, 58, 0, 23, 22, 22, 23, 59, 23, 60, 59, 0, 24, 23, 23, 24, 60, 24, 61, 60, 0, 25, 24, 24, 25, 61, 25, 62, 61, 0, 26, 25, 25, 26, 62, 26, 63, 62, 0, 27, 26, 26, 27, 63, 27, 64, 63, 0, 28, 27, 27, 28, 64, 28, 65, 64, 0, 29, 28, 28, 29, 65, 29, 66, 65, 0, 30, 29, 29, 30, 66, 30, 67, 66, 0, 31, 30, 30, 31, 67, 31, 68, 67, 0, 32, 31, 31, 32, 68, 32, 69, 68, 0, 33, 32, 32, 33, 69, 33, 70, 69, 0, 34, 33, 33, 34, 70, 34, 71, 70, 0, 35, 34, 34, 35, 71, 35, 72, 71, 0, 36, 35, 35, 36, 72, 36, 73, 72, 0, 1, 36, 36, 1, 73, 1, 38, 73, 37, 38, 39, 37, 39, 40, 37, 40, 41, 37, 41, 42, 37, 42, 43, 37, 43, 44, 37, 44, 45, 37, 45, 46, 37, 46, 47, 37, 47, 48, 37, 48, 49, 37, 49, 50, 37, 50, 51, 37, 51, 52, 37, 52, 53, 37, 53, 54, 37, 54, 55, 37, 55, 56, 37, 56, 57, 37, 57, 58, 37, 58, 59, 37, 59, 60, 37, 60, 61, 37, 61, 62, 37, 62, 63, 37, 63, 64, 37, 64, 65, 37, 65, 66, 37, 66, 67, 37, 67, 68, 37, 68, 69, 37, 69, 70, 37, 70, 71, 37, 71, 72, 37, 72, 73, 37, 73, 38], [ 0, 2, 1, 1, 2, 38, 2, 39, 38, 0, 3, 2, 2, 3, 39, 3, 40, 39, 0, 4, 3, 3, 4, 40, 4, 41, 40, 0, 5, 4, 4, 5, 41, 5, 42, 41, 0, 6, 5, 5, 6, 42, 6, 43, 42, 0, 7, 6, 6, 7, 43, 7, 44, 43, 0, 8, 7, 7, 8, 44, 8, 45, 44, 0, 9, 8, 8, 9, 45, 9, 46, 45, 0, 10, 9, 9, 10, 46, 10, 47, 46, 0, 11, 10, 10, 11, 47, 11, 48, 47, 0, 12, 11, 11, 12, 48, 12, 49, 48, 0, 13, 12, 12, 13, 49, 13, 50, 49, 0, 14, 13, 13, 14, 50, 14, 51, 50, 0, 15, 14, 14, 15, 51, 15, 52, 51, 0, 16, 15, 15, 16, 52, 16, 53, 52, 0, 17, 16, 16, 17, 53, 17, 54, 53, 0, 18, 17, 17, 18, 54, 18, 55, 54, 0, 19, 18, 18, 19, 55, 19, 56, 55, 0, 20, 19, 19, 20, 56, 20, 57, 56, 0, 21, 20, 20, 21, 57, 21, 58, 57, 0, 22, 21, 21, 22, 58, 22, 59, 58, 0, 23, 22, 22, 23, 59, 23, 60, 59, 0, 24, 23, 23, 24, 60, 24, 61, 60, 0, 25, 24, 24, 25, 61, 25, 62, 61, 0, 26, 25, 25, 26, 62, 26, 63, 62, 0, 27, 26, 26, 27, 63, 27, 64, 63, 0, 28, 27, 27, 28, 64, 28, 65, 64, 0, 29, 28, 28, 29, 65, 29, 66, 65, 0, 30, 29, 29, 30, 66, 30, 67, 66, 0, 31, 30, 30, 31, 67, 31, 68, 67, 0, 32, 31, 31, 32, 68, 32, 69, 68, 0, 33, 32, 32, 33, 69, 33, 70, 69, 0, 34, 33, 33, 34, 70, 34, 71, 70, 0, 35, 34, 34, 35, 71, 35, 72, 71, 0, 36, 35, 35, 36, 72, 36, 73, 72, 0, 1, 36, 36, 1, 73, 1, 38, 73, 37, 38, 39, 37, 39, 40, 37, 40, 41, 37, 41, 42, 37, 42, 43, 37, 43, 44, 37, 44, 45, 37, 45, 46, 37, 46, 47, 37, 47, 48, 37, 48, 49, 37, 49, 50, 37, 50, 51, 37, 51, 52, 37, 52, 53, 37, 53, 54, 37, 54, 55, 37, 55, 56, 37, 56, 57, 37, 57, 58, 37, 58, 59, 37, 59, 60, 37, 60, 61, 37, 61, 62, 37, 62, 63, 37, 63, 64, 37, 64, 65, 37, 65, 66, 37, 66, 67, 37, 67, 68, 37, 68, 69, 37, 69, 70, 37, 70, 71, 37, 71, 72, 37, 72, 73, 37, 73, 38], [ 0, 2, 1, 1, 2, 38, 2, 39, 38, 0, 3, 2, 2, 3, 39, 3, 40, 39, 0, 4, 3, 3, 4, 40, 4, 41, 40, 0, 5, 4, 4, 5, 41, 5, 42, 41, 0, 6, 5, 5, 6, 42, 6, 43, 42, 0, 7, 6, 6, 7, 43, 7, 44, 43, 0, 8, 7, 7, 8, 44, 8, 45, 44, 0, 9, 8, 8, 9, 45, 9, 46, 45, 0, 10, 9, 9, 10, 46, 10, 47, 46, 0, 11, 10, 10, 11, 47, 11, 48, 47, 0, 12, 11, 11, 12, 48, 12, 49, 48, 0, 13, 12, 12, 13, 49, 13, 50, 49, 0, 14, 13, 13, 14, 50, 14, 51, 50, 0, 15, 14, 14, 15, 51, 15, 52, 51, 0, 16, 15, 15, 16, 52, 16, 53, 52, 0, 17, 16, 16, 17, 53, 17, 54, 53, 0, 18, 17, 17, 18, 54, 18, 55, 54, 0, 19, 18, 18, 19, 55, 19, 56, 55, 0, 20, 19, 19, 20, 56, 20, 57, 56, 0, 21, 20, 20, 21, 57, 21, 58, 57, 0, 22, 21, 21, 22, 58, 22, 59, 58, 0, 23, 22, 22, 23, 59, 23, 60, 59, 0, 24, 23, 23, 24, 60, 24, 61, 60, 0, 25, 24, 24, 25, 61, 25, 62, 61, 0, 26, 25, 25, 26, 62, 26, 63, 62, 0, 27, 26, 26, 27, 63, 27, 64, 63, 0, 28, 27, 27, 28, 64, 28, 65, 64, 0, 29, 28, 28, 29, 65, 29, 66, 65, 0, 30, 29, 29, 30, 66, 30, 67, 66, 0, 31, 30, 30, 31, 67, 31, 68, 67, 0, 32, 31, 31, 32, 68, 32, 69, 68, 0, 33, 32, 32, 33, 69, 33, 70, 69, 0, 34, 33, 33, 34, 70, 34, 71, 70, 0, 35, 34, 34, 35, 71, 35, 72, 71, 0, 36, 35, 35, 36, 72, 36, 73, 72, 0, 1, 36, 36, 1, 73, 1, 38, 73, 37, 38, 39, 37, 39, 40, 37, 40, 41, 37, 41, 42, 37, 42, 43, 37, 43, 44, 37, 44, 45, 37, 45, 46, 37, 46, 47, 37, 47, 48, 37, 48, 49, 37, 49, 50, 37, 50, 51, 37, 51, 52, 37, 52, 53, 37, 53, 54, 37, 54, 55, 37, 55, 56, 37, 56, 57, 37, 57, 58, 37, 58, 59, 37, 59, 60, 37, 60, 61, 37, 61, 62, 37, 62, 63, 37, 63, 64, 37, 64, 65, 37, 65, 66, 37, 66, 67, 37, 67, 68, 37, 68, 69, 37, 69, 70, 37, 70, 71, 37, 71, 72, 37, 72, 73, 37, 73, 38], [ 0, 2, 1, 1, 2, 38, 2, 39, 38, 0, 3, 2, 2, 3, 39, 3, 40, 39, 0, 4, 3, 3, 4, 40, 4, 41, 40, 0, 5, 4, 4, 5, 41, 5, 42, 41, 0, 6, 5, 5, 6, 42, 6, 43, 42, 0, 7, 6, 6, 7, 43, 7, 44, 43, 0, 8, 7, 7, 8, 44, 8, 45, 44, 0, 9, 8, 8, 9, 45, 9, 46, 45, 0, 10, 9, 9, 10, 46, 10, 47, 46, 0, 11, 10, 10, 11, 47, 11, 48, 47, 0, 12, 11, 11, 12, 48, 12, 49, 48, 0, 13, 12, 12, 13, 49, 13, 50, 49, 0, 14, 13, 13, 14, 50, 14, 51, 50, 0, 15, 14, 14, 15, 51, 15, 52, 51, 0, 16, 15, 15, 16, 52, 16, 53, 52, 0, 17, 16, 16, 17, 53, 17, 54, 53, 0, 18, 17, 17, 18, 54, 18, 55, 54, 0, 19, 18, 18, 19, 55, 19, 56, 55, 0, 20, 19, 19, 20, 56, 20, 57, 56, 0, 21, 20, 20, 21, 57, 21, 58, 57, 0, 22, 21, 21, 22, 58, 22, 59, 58, 0, 23, 22, 22, 23, 59, 23, 60, 59, 0, 24, 23, 23, 24, 60, 24, 61, 60, 0, 25, 24, 24, 25, 61, 25, 62, 61, 0, 26, 25, 25, 26, 62, 26, 63, 62, 0, 27, 26, 26, 27, 63, 27, 64, 63, 0, 28, 27, 27, 28, 64, 28, 65, 64, 0, 29, 28, 28, 29, 65, 29, 66, 65, 0, 30, 29, 29, 30, 66, 30, 67, 66, 0, 31, 30, 30, 31, 67, 31, 68, 67, 0, 32, 31, 31, 32, 68, 32, 69, 68, 0, 33, 32, 32, 33, 69, 33, 70, 69, 0, 34, 33, 33, 34, 70, 34, 71, 70, 0, 35, 34, 34, 35, 71, 35, 72, 71, 0, 36, 35, 35, 36, 72, 36, 73, 72, 0, 1, 36, 36, 1, 73, 1, 38, 73, 37, 38, 39, 37, 39, 40, 37, 40, 41, 37, 41, 42, 37, 42, 43, 37, 43, 44, 37, 44, 45, 37, 45, 46, 37, 46, 47, 37, 47, 48, 37, 48, 49, 37, 49, 50, 37, 50, 51, 37, 51, 52, 37, 52, 53, 37, 53, 54, 37, 54, 55, 37, 55, 56, 37, 56, 57, 37, 57, 58, 37, 58, 59, 37, 59, 60, 37, 60, 61, 37, 61, 62, 37, 62, 63, 37, 63, 64, 37, 64, 65, 37, 65, 66, 37, 66, 67, 37, 67, 68, 37, 68, 69, 37, 69, 70, 37, 70, 71, 37, 71, 72, 37, 72, 73, 37, 73, 38], [ 0, 2, 1, 1, 2, 38, 2, 39, 38, 0, 3, 2, 2, 3, 39, 3, 40, 39, 0, 4, 3, 3, 4, 40, 4, 41, 40, 0, 5, 4, 4, 5, 41, 5, 42, 41, 0, 6, 5, 5, 6, 42, 6, 43, 42, 0, 7, 6, 6, 7, 43, 7, 44, 43, 0, 8, 7, 7, 8, 44, 8, 45, 44, 0, 9, 8, 8, 9, 45, 9, 46, 45, 0, 10, 9, 9, 10, 46, 10, 47, 46, 0, 11, 10, 10, 11, 47, 11, 48, 47, 0, 12, 11, 11, 12, 48, 12, 49, 48, 0, 13, 12, 12, 13, 49, 13, 50, 49, 0, 14, 13, 13, 14, 50, 14, 51, 50, 0, 15, 14, 14, 15, 51, 15, 52, 51, 0, 16, 15, 15, 16, 52, 16, 53, 52, 0, 17, 16, 16, 17, 53, 17, 54, 53, 0, 18, 17, 17, 18, 54, 18, 55, 54, 0, 19, 18, 18, 19, 55, 19, 56, 55, 0, 20, 19, 19, 20, 56, 20, 57, 56, 0, 21, 20, 20, 21, 57, 21, 58, 57, 0, 22, 21, 21, 22, 58, 22, 59, 58, 0, 23, 22, 22, 23, 59, 23, 60, 59, 0, 24, 23, 23, 24, 60, 24, 61, 60, 0, 25, 24, 24, 25, 61, 25, 62, 61, 0, 26, 25, 25, 26, 62, 26, 63, 62, 0, 27, 26, 26, 27, 63, 27, 64, 63, 0, 28, 27, 27, 28, 64, 28, 65, 64, 0, 29, 28, 28, 29, 65, 29, 66, 65, 0, 30, 29, 29, 30, 66, 30, 67, 66, 0, 31, 30, 30, 31, 67, 31, 68, 67, 0, 32, 31, 31, 32, 68, 32, 69, 68, 0, 33, 32, 32, 33, 69, 33, 70, 69, 0, 34, 33, 33, 34, 70, 34, 71, 70, 0, 35, 34, 34, 35, 71, 35, 72, 71, 0, 36, 35, 35, 36, 72, 36, 73, 72, 0, 1, 36, 36, 1, 73, 1, 38, 73, 37, 38, 39, 37, 39, 40, 37, 40, 41, 37, 41, 42, 37, 42, 43, 37, 43, 44, 37, 44, 45, 37, 45, 46, 37, 46, 47, 37, 47, 48, 37, 48, 49, 37, 49, 50, 37, 50, 51, 37, 51, 52, 37, 52, 53, 37, 53, 54, 37, 54, 55, 37, 55, 56, 37, 56, 57, 37, 57, 58, 37, 58, 59, 37, 59, 60, 37, 60, 61, 37, 61, 62, 37, 62, 63, 37, 63, 64, 37, 64, 65, 37, 65, 66, 37, 66, 67, 37, 67, 68, 37, 68, 69, 37, 69, 70, 37, 70, 71, 37, 71, 72, 37, 72, 73, 37, 73, 38], [ 0, 2, 1, 1, 2, 38, 2, 39, 38, 0, 3, 2, 2, 3, 39, 3, 40, 39, 0, 4, 3, 3, 4, 40, 4, 41, 40, 0, 5, 4, 4, 5, 41, 5, 42, 41, 0, 6, 5, 5, 6, 42, 6, 43, 42, 0, 7, 6, 6, 7, 43, 7, 44, 43, 0, 8, 7, 7, 8, 44, 8, 45, 44, 0, 9, 8, 8, 9, 45, 9, 46, 45, 0, 10, 9, 9, 10, 46, 10, 47, 46, 0, 11, 10, 10, 11, 47, 11, 48, 47, 0, 12, 11, 11, 12, 48, 12, 49, 48, 0, 13, 12, 12, 13, 49, 13, 50, 49, 0, 14, 13, 13, 14, 50, 14, 51, 50, 0, 15, 14, 14, 15, 51, 15, 52, 51, 0, 16, 15, 15, 16, 52, 16, 53, 52, 0, 17, 16, 16, 17, 53, 17, 54, 53, 0, 18, 17, 17, 18, 54, 18, 55, 54, 0, 19, 18, 18, 19, 55, 19, 56, 55, 0, 20, 19, 19, 20, 56, 20, 57, 56, 0, 21, 20, 20, 21, 57, 21, 58, 57, 0, 22, 21, 21, 22, 58, 22, 59, 58, 0, 23, 22, 22, 23, 59, 23, 60, 59, 0, 24, 23, 23, 24, 60, 24, 61, 60, 0, 25, 24, 24, 25, 61, 25, 62, 61, 0, 26, 25, 25, 26, 62, 26, 63, 62, 0, 27, 26, 26, 27, 63, 27, 64, 63, 0, 28, 27, 27, 28, 64, 28, 65, 64, 0, 29, 28, 28, 29, 65, 29, 66, 65, 0, 30, 29, 29, 30, 66, 30, 67, 66, 0, 31, 30, 30, 31, 67, 31, 68, 67, 0, 32, 31, 31, 32, 68, 32, 69, 68, 0, 33, 32, 32, 33, 69, 33, 70, 69, 0, 34, 33, 33, 34, 70, 34, 71, 70, 0, 35, 34, 34, 35, 71, 35, 72, 71, 0, 36, 35, 35, 36, 72, 36, 73, 72, 0, 1, 36, 36, 1, 73, 1, 38, 73, 37, 38, 39, 37, 39, 40, 37, 40, 41, 37, 41, 42, 37, 42, 43, 37, 43, 44, 37, 44, 45, 37, 45, 46, 37, 46, 47, 37, 47, 48, 37, 48, 49, 37, 49, 50, 37, 50, 51, 37, 51, 52, 37, 52, 53, 37, 53, 54, 37, 54, 55, 37, 55, 56, 37, 56, 57, 37, 57, 58, 37, 58, 59, 37, 59, 60, 37, 60, 61, 37, 61, 62, 37, 62, 63, 37, 63, 64, 37, 64, 65, 37, 65, 66, 37, 66, 67, 37, 67, 68, 37, 68, 69, 37, 69, 70, 37, 70, 71, 37, 71, 72, 37, 72, 73, 37, 73, 38], [ 0, 2, 1, 1, 2, 38, 2, 39, 38, 0, 3, 2, 2, 3, 39, 3, 40, 39, 0, 4, 3, 3, 4, 40, 4, 41, 40, 0, 5, 4, 4, 5, 41, 5, 42, 41, 0, 6, 5, 5, 6, 42, 6, 43, 42, 0, 7, 6, 6, 7, 43, 7, 44, 43, 0, 8, 7, 7, 8, 44, 8, 45, 44, 0, 9, 8, 8, 9, 45, 9, 46, 45, 0, 10, 9, 9, 10, 46, 10, 47, 46, 0, 11, 10, 10, 11, 47, 11, 48, 47, 0, 12, 11, 11, 12, 48, 12, 49, 48, 0, 13, 12, 12, 13, 49, 13, 50, 49, 0, 14, 13, 13, 14, 50, 14, 51, 50, 0, 15, 14, 14, 15, 51, 15, 52, 51, 0, 16, 15, 15, 16, 52, 16, 53, 52, 0, 17, 16, 16, 17, 53, 17, 54, 53, 0, 18, 17, 17, 18, 54, 18, 55, 54, 0, 19, 18, 18, 19, 55, 19, 56, 55, 0, 20, 19, 19, 20, 56, 20, 57, 56, 0, 21, 20, 20, 21, 57, 21, 58, 57, 0, 22, 21, 21, 22, 58, 22, 59, 58, 0, 23, 22, 22, 23, 59, 23, 60, 59, 0, 24, 23, 23, 24, 60, 24, 61, 60, 0, 25, 24, 24, 25, 61, 25, 62, 61, 0, 26, 25, 25, 26, 62, 26, 63, 62, 0, 27, 26, 26, 27, 63, 27, 64, 63, 0, 28, 27, 27, 28, 64, 28, 65, 64, 0, 29, 28, 28, 29, 65, 29, 66, 65, 0, 30, 29, 29, 30, 66, 30, 67, 66, 0, 31, 30, 30, 31, 67, 31, 68, 67, 0, 32, 31, 31, 32, 68, 32, 69, 68, 0, 33, 32, 32, 33, 69, 33, 70, 69, 0, 34, 33, 33, 34, 70, 34, 71, 70, 0, 35, 34, 34, 35, 71, 35, 72, 71, 0, 36, 35, 35, 36, 72, 36, 73, 72, 0, 1, 36, 36, 1, 73, 1, 38, 73, 37, 38, 39, 37, 39, 40, 37, 40, 41, 37, 41, 42, 37, 42, 43, 37, 43, 44, 37, 44, 45, 37, 45, 46, 37, 46, 47, 37, 47, 48, 37, 48, 49, 37, 49, 50, 37, 50, 51, 37, 51, 52, 37, 52, 53, 37, 53, 54, 37, 54, 55, 37, 55, 56, 37, 56, 57, 37, 57, 58, 37, 58, 59, 37, 59, 60, 37, 60, 61, 37, 61, 62, 37, 62, 63, 37, 63, 64, 37, 64, 65, 37, 65, 66, 37, 66, 67, 37, 67, 68, 37, 68, 69, 37, 69, 70, 37, 70, 71, 37, 71, 72, 37, 72, 73, 37, 73, 38], [ 0, 2, 1, 1, 2, 38, 2, 39, 38, 0, 3, 2, 2, 3, 39, 3, 40, 39, 0, 4, 3, 3, 4, 40, 4, 41, 40, 0, 5, 4, 4, 5, 41, 5, 42, 41, 0, 6, 5, 5, 6, 42, 6, 43, 42, 0, 7, 6, 6, 7, 43, 7, 44, 43, 0, 8, 7, 7, 8, 44, 8, 45, 44, 0, 9, 8, 8, 9, 45, 9, 46, 45, 0, 10, 9, 9, 10, 46, 10, 47, 46, 0, 11, 10, 10, 11, 47, 11, 48, 47, 0, 12, 11, 11, 12, 48, 12, 49, 48, 0, 13, 12, 12, 13, 49, 13, 50, 49, 0, 14, 13, 13, 14, 50, 14, 51, 50, 0, 15, 14, 14, 15, 51, 15, 52, 51, 0, 16, 15, 15, 16, 52, 16, 53, 52, 0, 17, 16, 16, 17, 53, 17, 54, 53, 0, 18, 17, 17, 18, 54, 18, 55, 54, 0, 19, 18, 18, 19, 55, 19, 56, 55, 0, 20, 19, 19, 20, 56, 20, 57, 56, 0, 21, 20, 20, 21, 57, 21, 58, 57, 0, 22, 21, 21, 22, 58, 22, 59, 58, 0, 23, 22, 22, 23, 59, 23, 60, 59, 0, 24, 23, 23, 24, 60, 24, 61, 60, 0, 25, 24, 24, 25, 61, 25, 62, 61, 0, 26, 25, 25, 26, 62, 26, 63, 62, 0, 27, 26, 26, 27, 63, 27, 64, 63, 0, 28, 27, 27, 28, 64, 28, 65, 64, 0, 29, 28, 28, 29, 65, 29, 66, 65, 0, 30, 29, 29, 30, 66, 30, 67, 66, 0, 31, 30, 30, 31, 67, 31, 68, 67, 0, 32, 31, 31, 32, 68, 32, 69, 68, 0, 33, 32, 32, 33, 69, 33, 70, 69, 0, 34, 33, 33, 34, 70, 34, 71, 70, 0, 35, 34, 34, 35, 71, 35, 72, 71, 0, 36, 35, 35, 36, 72, 36, 73, 72, 0, 1, 36, 36, 1, 73, 1, 38, 73, 37, 38, 39, 37, 39, 40, 37, 40, 41, 37, 41, 42, 37, 42, 43, 37, 43, 44, 37, 44, 45, 37, 45, 46, 37, 46, 47, 37, 47, 48, 37, 48, 49, 37, 49, 50, 37, 50, 51, 37, 51, 52, 37, 52, 53, 37, 53, 54, 37, 54, 55, 37, 55, 56, 37, 56, 57, 37, 57, 58, 37, 58, 59, 37, 59, 60, 37, 60, 61, 37, 61, 62, 37, 62, 63, 37, 63, 64, 37, 64, 65, 37, 65, 66, 37, 66, 67, 37, 67, 68, 37, 68, 69, 37, 69, 70, 37, 70, 71, 37, 71, 72, 37, 72, 73, 37, 73, 38], [ 0, 2, 1, 1, 2, 38, 2, 39, 38, 0, 3, 2, 2, 3, 39, 3, 40, 39, 0, 4, 3, 3, 4, 40, 4, 41, 40, 0, 5, 4, 4, 5, 41, 5, 42, 41, 0, 6, 5, 5, 6, 42, 6, 43, 42, 0, 7, 6, 6, 7, 43, 7, 44, 43, 0, 8, 7, 7, 8, 44, 8, 45, 44, 0, 9, 8, 8, 9, 45, 9, 46, 45, 0, 10, 9, 9, 10, 46, 10, 47, 46, 0, 11, 10, 10, 11, 47, 11, 48, 47, 0, 12, 11, 11, 12, 48, 12, 49, 48, 0, 13, 12, 12, 13, 49, 13, 50, 49, 0, 14, 13, 13, 14, 50, 14, 51, 50, 0, 15, 14, 14, 15, 51, 15, 52, 51, 0, 16, 15, 15, 16, 52, 16, 53, 52, 0, 17, 16, 16, 17, 53, 17, 54, 53, 0, 18, 17, 17, 18, 54, 18, 55, 54, 0, 19, 18, 18, 19, 55, 19, 56, 55, 0, 20, 19, 19, 20, 56, 20, 57, 56, 0, 21, 20, 20, 21, 57, 21, 58, 57, 0, 22, 21, 21, 22, 58, 22, 59, 58, 0, 23, 22, 22, 23, 59, 23, 60, 59, 0, 24, 23, 23, 24, 60, 24, 61, 60, 0, 25, 24, 24, 25, 61, 25, 62, 61, 0, 26, 25, 25, 26, 62, 26, 63, 62, 0, 27, 26, 26, 27, 63, 27, 64, 63, 0, 28, 27, 27, 28, 64, 28, 65, 64, 0, 29, 28, 28, 29, 65, 29, 66, 65, 0, 30, 29, 29, 30, 66, 30, 67, 66, 0, 31, 30, 30, 31, 67, 31, 68, 67, 0, 32, 31, 31, 32, 68, 32, 69, 68, 0, 33, 32, 32, 33, 69, 33, 70, 69, 0, 34, 33, 33, 34, 70, 34, 71, 70, 0, 35, 34, 34, 35, 71, 35, 72, 71, 0, 36, 35, 35, 36, 72, 36, 73, 72, 0, 1, 36, 36, 1, 73, 1, 38, 73, 37, 38, 39, 37, 39, 40, 37, 40, 41, 37, 41, 42, 37, 42, 43, 37, 43, 44, 37, 44, 45, 37, 45, 46, 37, 46, 47, 37, 47, 48, 37, 48, 49, 37, 49, 50, 37, 50, 51, 37, 51, 52, 37, 52, 53, 37, 53, 54, 37, 54, 55, 37, 55, 56, 37, 56, 57, 37, 57, 58, 37, 58, 59, 37, 59, 60, 37, 60, 61, 37, 61, 62, 37, 62, 63, 37, 63, 64, 37, 64, 65, 37, 65, 66, 37, 66, 67, 37, 67, 68, 37, 68, 69, 37, 69, 70, 37, 70, 71, 37, 71, 72, 37, 72, 73, 37, 73, 38], [ 0, 2, 1, 1, 2, 38, 2, 39, 38, 0, 3, 2, 2, 3, 39, 3, 40, 39, 0, 4, 3, 3, 4, 40, 4, 41, 40, 0, 5, 4, 4, 5, 41, 5, 42, 41, 0, 6, 5, 5, 6, 42, 6, 43, 42, 0, 7, 6, 6, 7, 43, 7, 44, 43, 0, 8, 7, 7, 8, 44, 8, 45, 44, 0, 9, 8, 8, 9, 45, 9, 46, 45, 0, 10, 9, 9, 10, 46, 10, 47, 46, 0, 11, 10, 10, 11, 47, 11, 48, 47, 0, 12, 11, 11, 12, 48, 12, 49, 48, 0, 13, 12, 12, 13, 49, 13, 50, 49, 0, 14, 13, 13, 14, 50, 14, 51, 50, 0, 15, 14, 14, 15, 51, 15, 52, 51, 0, 16, 15, 15, 16, 52, 16, 53, 52, 0, 17, 16, 16, 17, 53, 17, 54, 53, 0, 18, 17, 17, 18, 54, 18, 55, 54, 0, 19, 18, 18, 19, 55, 19, 56, 55, 0, 20, 19, 19, 20, 56, 20, 57, 56, 0, 21, 20, 20, 21, 57, 21, 58, 57, 0, 22, 21, 21, 22, 58, 22, 59, 58, 0, 23, 22, 22, 23, 59, 23, 60, 59, 0, 24, 23, 23, 24, 60, 24, 61, 60, 0, 25, 24, 24, 25, 61, 25, 62, 61, 0, 26, 25, 25, 26, 62, 26, 63, 62, 0, 27, 26, 26, 27, 63, 27, 64, 63, 0, 28, 27, 27, 28, 64, 28, 65, 64, 0, 29, 28, 28, 29, 65, 29, 66, 65, 0, 30, 29, 29, 30, 66, 30, 67, 66, 0, 31, 30, 30, 31, 67, 31, 68, 67, 0, 32, 31, 31, 32, 68, 32, 69, 68, 0, 33, 32, 32, 33, 69, 33, 70, 69, 0, 34, 33, 33, 34, 70, 34, 71, 70, 0, 35, 34, 34, 35, 71, 35, 72, 71, 0, 36, 35, 35, 36, 72, 36, 73, 72, 0, 1, 36, 36, 1, 73, 1, 38, 73, 37, 38, 39, 37, 39, 40, 37, 40, 41, 37, 41, 42, 37, 42, 43, 37, 43, 44, 37, 44, 45, 37, 45, 46, 37, 46, 47, 37, 47, 48, 37, 48, 49, 37, 49, 50, 37, 50, 51, 37, 51, 52, 37, 52, 53, 37, 53, 54, 37, 54, 55, 37, 55, 56, 37, 56, 57, 37, 57, 58, 37, 58, 59, 37, 59, 60, 37, 60, 61, 37, 61, 62, 37, 62, 63, 37, 63, 64, 37, 64, 65, 37, 65, 66, 37, 66, 67, 37, 67, 68, 37, 68, 69, 37, 69, 70, 37, 70, 71, 37, 71, 72, 37, 72, 73, 37, 73, 38], [ 0, 2, 1, 1, 2, 38, 2, 39, 38, 0, 3, 2, 2, 3, 39, 3, 40, 39, 0, 4, 3, 3, 4, 40, 4, 41, 40, 0, 5, 4, 4, 5, 41, 5, 42, 41, 0, 6, 5, 5, 6, 42, 6, 43, 42, 0, 7, 6, 6, 7, 43, 7, 44, 43, 0, 8, 7, 7, 8, 44, 8, 45, 44, 0, 9, 8, 8, 9, 45, 9, 46, 45, 0, 10, 9, 9, 10, 46, 10, 47, 46, 0, 11, 10, 10, 11, 47, 11, 48, 47, 0, 12, 11, 11, 12, 48, 12, 49, 48, 0, 13, 12, 12, 13, 49, 13, 50, 49, 0, 14, 13, 13, 14, 50, 14, 51, 50, 0, 15, 14, 14, 15, 51, 15, 52, 51, 0, 16, 15, 15, 16, 52, 16, 53, 52, 0, 17, 16, 16, 17, 53, 17, 54, 53, 0, 18, 17, 17, 18, 54, 18, 55, 54, 0, 19, 18, 18, 19, 55, 19, 56, 55, 0, 20, 19, 19, 20, 56, 20, 57, 56, 0, 21, 20, 20, 21, 57, 21, 58, 57, 0, 22, 21, 21, 22, 58, 22, 59, 58, 0, 23, 22, 22, 23, 59, 23, 60, 59, 0, 24, 23, 23, 24, 60, 24, 61, 60, 0, 25, 24, 24, 25, 61, 25, 62, 61, 0, 26, 25, 25, 26, 62, 26, 63, 62, 0, 27, 26, 26, 27, 63, 27, 64, 63, 0, 28, 27, 27, 28, 64, 28, 65, 64, 0, 29, 28, 28, 29, 65, 29, 66, 65, 0, 30, 29, 29, 30, 66, 30, 67, 66, 0, 31, 30, 30, 31, 67, 31, 68, 67, 0, 32, 31, 31, 32, 68, 32, 69, 68, 0, 33, 32, 32, 33, 69, 33, 70, 69, 0, 34, 33, 33, 34, 70, 34, 71, 70, 0, 35, 34, 34, 35, 71, 35, 72, 71, 0, 36, 35, 35, 36, 72, 36, 73, 72, 0, 1, 36, 36, 1, 73, 1, 38, 73, 37, 38, 39, 37, 39, 40, 37, 40, 41, 37, 41, 42, 37, 42, 43, 37, 43, 44, 37, 44, 45, 37, 45, 46, 37, 46, 47, 37, 47, 48, 37, 48, 49, 37, 49, 50, 37, 50, 51, 37, 51, 52, 37, 52, 53, 37, 53, 54, 37, 54, 55, 37, 55, 56, 37, 56, 57, 37, 57, 58, 37, 58, 59, 37, 59, 60, 37, 60, 61, 37, 61, 62, 37, 62, 63, 37, 63, 64, 37, 64, 65, 37, 65, 66, 37, 66, 67, 37, 67, 68, 37, 68, 69, 37, 69, 70, 37, 70, 71, 37, 71, 72, 37, 72, 73, 37, 73, 38], [ 0, 2, 1, 1, 2, 38, 2, 39, 38, 0, 3, 2, 2, 3, 39, 3, 40, 39, 0, 4, 3, 3, 4, 40, 4, 41, 40, 0, 5, 4, 4, 5, 41, 5, 42, 41, 0, 6, 5, 5, 6, 42, 6, 43, 42, 0, 7, 6, 6, 7, 43, 7, 44, 43, 0, 8, 7, 7, 8, 44, 8, 45, 44, 0, 9, 8, 8, 9, 45, 9, 46, 45, 0, 10, 9, 9, 10, 46, 10, 47, 46, 0, 11, 10, 10, 11, 47, 11, 48, 47, 0, 12, 11, 11, 12, 48, 12, 49, 48, 0, 13, 12, 12, 13, 49, 13, 50, 49, 0, 14, 13, 13, 14, 50, 14, 51, 50, 0, 15, 14, 14, 15, 51, 15, 52, 51, 0, 16, 15, 15, 16, 52, 16, 53, 52, 0, 17, 16, 16, 17, 53, 17, 54, 53, 0, 18, 17, 17, 18, 54, 18, 55, 54, 0, 19, 18, 18, 19, 55, 19, 56, 55, 0, 20, 19, 19, 20, 56, 20, 57, 56, 0, 21, 20, 20, 21, 57, 21, 58, 57, 0, 22, 21, 21, 22, 58, 22, 59, 58, 0, 23, 22, 22, 23, 59, 23, 60, 59, 0, 24, 23, 23, 24, 60, 24, 61, 60, 0, 25, 24, 24, 25, 61, 25, 62, 61, 0, 26, 25, 25, 26, 62, 26, 63, 62, 0, 27, 26, 26, 27, 63, 27, 64, 63, 0, 28, 27, 27, 28, 64, 28, 65, 64, 0, 29, 28, 28, 29, 65, 29, 66, 65, 0, 30, 29, 29, 30, 66, 30, 67, 66, 0, 31, 30, 30, 31, 67, 31, 68, 67, 0, 32, 31, 31, 32, 68, 32, 69, 68, 0, 33, 32, 32, 33, 69, 33, 70, 69, 0, 34, 33, 33, 34, 70, 34, 71, 70, 0, 35, 34, 34, 35, 71, 35, 72, 71, 0, 36, 35, 35, 36, 72, 36, 73, 72, 0, 1, 36, 36, 1, 73, 1, 38, 73, 37, 38, 39, 37, 39, 40, 37, 40, 41, 37, 41, 42, 37, 42, 43, 37, 43, 44, 37, 44, 45, 37, 45, 46, 37, 46, 47, 37, 47, 48, 37, 48, 49, 37, 49, 50, 37, 50, 51, 37, 51, 52, 37, 52, 53, 37, 53, 54, 37, 54, 55, 37, 55, 56, 37, 56, 57, 37, 57, 58, 37, 58, 59, 37, 59, 60, 37, 60, 61, 37, 61, 62, 37, 62, 63, 37, 63, 64, 37, 64, 65, 37, 65, 66, 37, 66, 67, 37, 67, 68, 37, 68, 69, 37, 69, 70, 37, 70, 71, 37, 71, 72, 37, 72, 73, 37, 73, 38], [ 2, 1, 0, 2, 4, 3, 2, 5, 4, 2, 0, 37, 2, 37, 5, 37, 36, 5, 5, 36, 6, 21, 20, 19, 21, 23, 22, 21, 24, 23, 21, 19, 18, 21, 18, 24, 18, 17, 24, 24, 17, 25, 40, 38, 39, 40, 41, 42, 40, 42, 43, 40, 75, 38, 40, 43, 75, 75, 43, 74, 43, 44, 74, 59, 57, 58, 59, 60, 61, 59, 61, 62, 59, 56, 57, 59, 62, 56, 56, 62, 55, 62, 63, 55, 6, 36, 7, 36, 35, 7, 44, 45, 74, 74, 45, 73, 7, 35, 8, 35, 34, 8, 45, 46, 73, 73, 46, 72, 8, 34, 9, 34, 33, 9, 46, 47, 72, 72, 47, 71, 9, 33, 10, 33, 32, 10, 47, 48, 71, 71, 48, 70, 10, 32, 11, 32, 31, 11, 48, 49, 70, 70, 49, 69, 11, 31, 12, 31, 30, 12, 49, 50, 69, 69, 50, 68, 12, 30, 13, 30, 29, 13, 50, 51, 68, 68, 51, 67, 13, 29, 14, 29, 28, 14, 51, 52, 67, 67, 52, 66, 14, 28, 15, 28, 27, 15, 52, 53, 66, 66, 53, 65, 15, 27, 16, 27, 26, 16, 53, 54, 65, 65, 54, 64, 16, 26, 17, 26, 25, 17, 54, 55, 64, 64, 55, 63, 0, 1, 39, 0, 39, 38, 1, 2, 40, 1, 40, 39, 2, 3, 41, 2, 41, 40, 3, 4, 42, 3, 42, 41, 6, 7, 45, 6, 45, 44, 7, 8, 46, 7, 46, 45, 8, 9, 47, 8, 47, 46, 9, 10, 48, 9, 48, 47, 10, 11, 49, 10, 49, 48, 11, 12, 50, 11, 50, 49, 12, 13, 51, 12, 51, 50, 13, 14, 52, 13, 52, 51, 14, 15, 53, 14, 53, 52, 15, 16, 54, 15, 54, 53, 16, 17, 55, 16, 55, 54, 19, 20, 58, 19, 58, 57, 20, 21, 59, 20, 59, 58, 21, 22, 60, 21, 60, 59, 22, 23, 61, 22, 61, 60, 25, 26, 64, 25, 64, 63, 26, 27, 65, 26, 65, 64, 27, 28, 66, 27, 66, 65, 28, 29, 67, 28, 67, 66, 29, 30, 68, 29, 68, 67, 30, 31, 69, 30, 69, 68, 31, 32, 70, 31, 70, 69, 32, 33, 71, 32, 71, 70, 33, 34, 72, 33, 72, 71, 34, 35, 73, 34, 73, 72, 35, 36, 74, 35, 74, 73, 37, 0, 36, 76, 36, 0, 76, 74, 36, 76, 38, 74, 75, 74, 38, 4, 5, 6, 4, 6, 77, 77, 6, 44, 42, 77, 44, 42, 44, 43, 18, 19, 17, 17, 19, 78, 17, 78, 55, 55, 78, 57, 56, 55, 57, 25, 23, 24, 25, 79, 23, 79, 25, 63, 61, 79, 63, 61, 63, 62], [ 0, 2, 1, 1, 2, 38, 2, 39, 38, 0, 3, 2, 2, 3, 39, 3, 40, 39, 0, 4, 3, 3, 4, 40, 4, 41, 40, 0, 5, 4, 4, 5, 41, 5, 42, 41, 0, 6, 5, 5, 6, 42, 6, 43, 42, 0, 7, 6, 6, 7, 43, 7, 44, 43, 0, 8, 7, 7, 8, 44, 8, 45, 44, 0, 9, 8, 8, 9, 45, 9, 46, 45, 0, 10, 9, 9, 10, 46, 10, 47, 46, 0, 11, 10, 10, 11, 47, 11, 48, 47, 0, 12, 11, 11, 12, 48, 12, 49, 48, 0, 13, 12, 12, 13, 49, 13, 50, 49, 0, 14, 13, 13, 14, 50, 14, 51, 50, 0, 15, 14, 14, 15, 51, 15, 52, 51, 0, 16, 15, 15, 16, 52, 16, 53, 52, 0, 17, 16, 16, 17, 53, 17, 54, 53, 0, 18, 17, 17, 18, 54, 18, 55, 54, 0, 19, 18, 18, 19, 55, 19, 56, 55, 0, 20, 19, 19, 20, 56, 20, 57, 56, 0, 21, 20, 20, 21, 57, 21, 58, 57, 0, 22, 21, 21, 22, 58, 22, 59, 58, 0, 23, 22, 22, 23, 59, 23, 60, 59, 0, 24, 23, 23, 24, 60, 24, 61, 60, 0, 25, 24, 24, 25, 61, 25, 62, 61, 0, 26, 25, 25, 26, 62, 26, 63, 62, 0, 27, 26, 26, 27, 63, 27, 64, 63, 0, 28, 27, 27, 28, 64, 28, 65, 64, 0, 29, 28, 28, 29, 65, 29, 66, 65, 0, 30, 29, 29, 30, 66, 30, 67, 66, 0, 31, 30, 30, 31, 67, 31, 68, 67, 0, 32, 31, 31, 32, 68, 32, 69, 68, 0, 33, 32, 32, 33, 69, 33, 70, 69, 0, 34, 33, 33, 34, 70, 34, 71, 70, 0, 35, 34, 34, 35, 71, 35, 72, 71, 0, 36, 35, 35, 36, 72, 36, 73, 72, 0, 1, 36, 36, 1, 73, 1, 38, 73, 37, 38, 39, 37, 39, 40, 37, 40, 41, 37, 41, 42, 37, 42, 43, 37, 43, 44, 37, 44, 45, 37, 45, 46, 37, 46, 47, 37, 47, 48, 37, 48, 49, 37, 49, 50, 37, 50, 51, 37, 51, 52, 37, 52, 53, 37, 53, 54, 37, 54, 55, 37, 55, 56, 37, 56, 57, 37, 57, 58, 37, 58, 59, 37, 59, 60, 37, 60, 61, 37, 61, 62, 37, 62, 63, 37, 63, 64, 37, 64, 65, 37, 65, 66, 37, 66, 67, 37, 67, 68, 37, 68, 69, 37, 69, 70, 37, 70, 71, 37, 71, 72, 37, 72, 73, 37, 73, 38], [ 0, 2, 1, 1, 2, 38, 2, 39, 38, 0, 3, 2, 2, 3, 39, 3, 40, 39, 0, 4, 3, 3, 4, 40, 4, 41, 40, 0, 5, 4, 4, 5, 41, 5, 42, 41, 0, 6, 5, 5, 6, 42, 6, 43, 42, 0, 7, 6, 6, 7, 43, 7, 44, 43, 0, 8, 7, 7, 8, 44, 8, 45, 44, 0, 9, 8, 8, 9, 45, 9, 46, 45, 0, 10, 9, 9, 10, 46, 10, 47, 46, 0, 11, 10, 10, 11, 47, 11, 48, 47, 0, 12, 11, 11, 12, 48, 12, 49, 48, 0, 13, 12, 12, 13, 49, 13, 50, 49, 0, 14, 13, 13, 14, 50, 14, 51, 50, 0, 15, 14, 14, 15, 51, 15, 52, 51, 0, 16, 15, 15, 16, 52, 16, 53, 52, 0, 17, 16, 16, 17, 53, 17, 54, 53, 0, 18, 17, 17, 18, 54, 18, 55, 54, 0, 19, 18, 18, 19, 55, 19, 56, 55, 0, 20, 19, 19, 20, 56, 20, 57, 56, 0, 21, 20, 20, 21, 57, 21, 58, 57, 0, 22, 21, 21, 22, 58, 22, 59, 58, 0, 23, 22, 22, 23, 59, 23, 60, 59, 0, 24, 23, 23, 24, 60, 24, 61, 60, 0, 25, 24, 24, 25, 61, 25, 62, 61, 0, 26, 25, 25, 26, 62, 26, 63, 62, 0, 27, 26, 26, 27, 63, 27, 64, 63, 0, 28, 27, 27, 28, 64, 28, 65, 64, 0, 29, 28, 28, 29, 65, 29, 66, 65, 0, 30, 29, 29, 30, 66, 30, 67, 66, 0, 31, 30, 30, 31, 67, 31, 68, 67, 0, 32, 31, 31, 32, 68, 32, 69, 68, 0, 33, 32, 32, 33, 69, 33, 70, 69, 0, 34, 33, 33, 34, 70, 34, 71, 70, 0, 35, 34, 34, 35, 71, 35, 72, 71, 0, 36, 35, 35, 36, 72, 36, 73, 72, 0, 1, 36, 36, 1, 73, 1, 38, 73, 37, 38, 39, 37, 39, 40, 37, 40, 41, 37, 41, 42, 37, 42, 43, 37, 43, 44, 37, 44, 45, 37, 45, 46, 37, 46, 47, 37, 47, 48, 37, 48, 49, 37, 49, 50, 37, 50, 51, 37, 51, 52, 37, 52, 53, 37, 53, 54, 37, 54, 55, 37, 55, 56, 37, 56, 57, 37, 57, 58, 37, 58, 59, 37, 59, 60, 37, 60, 61, 37, 61, 62, 37, 62, 63, 37, 63, 64, 37, 64, 65, 37, 65, 66, 37, 66, 67, 37, 67, 68, 37, 68, 69, 37, 69, 70, 37, 70, 71, 37, 71, 72, 37, 72, 73, 37, 73, 38], [ 0, 2, 1, 1, 2, 38, 2, 39, 38, 0, 3, 2, 2, 3, 39, 3, 40, 39, 0, 4, 3, 3, 4, 40, 4, 41, 40, 0, 5, 4, 4, 5, 41, 5, 42, 41, 0, 6, 5, 5, 6, 42, 6, 43, 42, 0, 7, 6, 6, 7, 43, 7, 44, 43, 0, 8, 7, 7, 8, 44, 8, 45, 44, 0, 9, 8, 8, 9, 45, 9, 46, 45, 0, 10, 9, 9, 10, 46, 10, 47, 46, 0, 11, 10, 10, 11, 47, 11, 48, 47, 0, 12, 11, 11, 12, 48, 12, 49, 48, 0, 13, 12, 12, 13, 49, 13, 50, 49, 0, 14, 13, 13, 14, 50, 14, 51, 50, 0, 15, 14, 14, 15, 51, 15, 52, 51, 0, 16, 15, 15, 16, 52, 16, 53, 52, 0, 17, 16, 16, 17, 53, 17, 54, 53, 0, 18, 17, 17, 18, 54, 18, 55, 54, 0, 19, 18, 18, 19, 55, 19, 56, 55, 0, 20, 19, 19, 20, 56, 20, 57, 56, 0, 21, 20, 20, 21, 57, 21, 58, 57, 0, 22, 21, 21, 22, 58, 22, 59, 58, 0, 23, 22, 22, 23, 59, 23, 60, 59, 0, 24, 23, 23, 24, 60, 24, 61, 60, 0, 25, 24, 24, 25, 61, 25, 62, 61, 0, 26, 25, 25, 26, 62, 26, 63, 62, 0, 27, 26, 26, 27, 63, 27, 64, 63, 0, 28, 27, 27, 28, 64, 28, 65, 64, 0, 29, 28, 28, 29, 65, 29, 66, 65, 0, 30, 29, 29, 30, 66, 30, 67, 66, 0, 31, 30, 30, 31, 67, 31, 68, 67, 0, 32, 31, 31, 32, 68, 32, 69, 68, 0, 33, 32, 32, 33, 69, 33, 70, 69, 0, 34, 33, 33, 34, 70, 34, 71, 70, 0, 35, 34, 34, 35, 71, 35, 72, 71, 0, 36, 35, 35, 36, 72, 36, 73, 72, 0, 1, 36, 36, 1, 73, 1, 38, 73, 37, 38, 39, 37, 39, 40, 37, 40, 41, 37, 41, 42, 37, 42, 43, 37, 43, 44, 37, 44, 45, 37, 45, 46, 37, 46, 47, 37, 47, 48, 37, 48, 49, 37, 49, 50, 37, 50, 51, 37, 51, 52, 37, 52, 53, 37, 53, 54, 37, 54, 55, 37, 55, 56, 37, 56, 57, 37, 57, 58, 37, 58, 59, 37, 59, 60, 37, 60, 61, 37, 61, 62, 37, 62, 63, 37, 63, 64, 37, 64, 65, 37, 65, 66, 37, 66, 67, 37, 67, 68, 37, 68, 69, 37, 69, 70, 37, 70, 71, 37, 71, 72, 37, 72, 73, 37, 73, 38], [ 0, 2, 1, 1, 2, 38, 2, 39, 38, 0, 3, 2, 2, 3, 39, 3, 40, 39, 0, 4, 3, 3, 4, 40, 4, 41, 40, 0, 5, 4, 4, 5, 41, 5, 42, 41, 0, 6, 5, 5, 6, 42, 6, 43, 42, 0, 7, 6, 6, 7, 43, 7, 44, 43, 0, 8, 7, 7, 8, 44, 8, 45, 44, 0, 9, 8, 8, 9, 45, 9, 46, 45, 0, 10, 9, 9, 10, 46, 10, 47, 46, 0, 11, 10, 10, 11, 47, 11, 48, 47, 0, 12, 11, 11, 12, 48, 12, 49, 48, 0, 13, 12, 12, 13, 49, 13, 50, 49, 0, 14, 13, 13, 14, 50, 14, 51, 50, 0, 15, 14, 14, 15, 51, 15, 52, 51, 0, 16, 15, 15, 16, 52, 16, 53, 52, 0, 17, 16, 16, 17, 53, 17, 54, 53, 0, 18, 17, 17, 18, 54, 18, 55, 54, 0, 19, 18, 18, 19, 55, 19, 56, 55, 0, 20, 19, 19, 20, 56, 20, 57, 56, 0, 21, 20, 20, 21, 57, 21, 58, 57, 0, 22, 21, 21, 22, 58, 22, 59, 58, 0, 23, 22, 22, 23, 59, 23, 60, 59, 0, 24, 23, 23, 24, 60, 24, 61, 60, 0, 25, 24, 24, 25, 61, 25, 62, 61, 0, 26, 25, 25, 26, 62, 26, 63, 62, 0, 27, 26, 26, 27, 63, 27, 64, 63, 0, 28, 27, 27, 28, 64, 28, 65, 64, 0, 29, 28, 28, 29, 65, 29, 66, 65, 0, 30, 29, 29, 30, 66, 30, 67, 66, 0, 31, 30, 30, 31, 67, 31, 68, 67, 0, 32, 31, 31, 32, 68, 32, 69, 68, 0, 33, 32, 32, 33, 69, 33, 70, 69, 0, 34, 33, 33, 34, 70, 34, 71, 70, 0, 35, 34, 34, 35, 71, 35, 72, 71, 0, 36, 35, 35, 36, 72, 36, 73, 72, 0, 1, 36, 36, 1, 73, 1, 38, 73, 37, 38, 39, 37, 39, 40, 37, 40, 41, 37, 41, 42, 37, 42, 43, 37, 43, 44, 37, 44, 45, 37, 45, 46, 37, 46, 47, 37, 47, 48, 37, 48, 49, 37, 49, 50, 37, 50, 51, 37, 51, 52, 37, 52, 53, 37, 53, 54, 37, 54, 55, 37, 55, 56, 37, 56, 57, 37, 57, 58, 37, 58, 59, 37, 59, 60, 37, 60, 61, 37, 61, 62, 37, 62, 63, 37, 63, 64, 37, 64, 65, 37, 65, 66, 37, 66, 67, 37, 67, 68, 37, 68, 69, 37, 69, 70, 37, 70, 71, 37, 71, 72, 37, 72, 73, 37, 73, 38], [ 2, 1, 0, 2, 4, 3, 2, 5, 4, 2, 0, 37, 2, 37, 5, 37, 36, 5, 5, 36, 6, 21, 20, 19, 21, 23, 22, 21, 24, 23, 21, 19, 18, 21, 18, 24, 18, 17, 24, 24, 17, 25, 40, 38, 39, 40, 41, 42, 40, 42, 43, 40, 75, 38, 40, 43, 75, 75, 43, 74, 43, 44, 74, 59, 57, 58, 59, 60, 61, 59, 61, 62, 59, 56, 57, 59, 62, 56, 56, 62, 55, 62, 63, 55, 6, 36, 7, 36, 35, 7, 44, 45, 74, 74, 45, 73, 7, 35, 8, 35, 34, 8, 45, 46, 73, 73, 46, 72, 8, 34, 9, 34, 33, 9, 46, 47, 72, 72, 47, 71, 9, 33, 10, 33, 32, 10, 47, 48, 71, 71, 48, 70, 10, 32, 11, 32, 31, 11, 48, 49, 70, 70, 49, 69, 11, 31, 12, 31, 30, 12, 49, 50, 69, 69, 50, 68, 12, 30, 13, 30, 29, 13, 50, 51, 68, 68, 51, 67, 13, 29, 14, 29, 28, 14, 51, 52, 67, 67, 52, 66, 14, 28, 15, 28, 27, 15, 52, 53, 66, 66, 53, 65, 15, 27, 16, 27, 26, 16, 53, 54, 65, 65, 54, 64, 16, 26, 17, 26, 25, 17, 54, 55, 64, 64, 55, 63, 0, 1, 39, 0, 39, 38, 1, 2, 40, 1, 40, 39, 2, 3, 41, 2, 41, 40, 3, 4, 42, 3, 42, 41, 6, 7, 45, 6, 45, 44, 7, 8, 46, 7, 46, 45, 8, 9, 47, 8, 47, 46, 9, 10, 48, 9, 48, 47, 10, 11, 49, 10, 49, 48, 11, 12, 50, 11, 50, 49, 12, 13, 51, 12, 51, 50, 13, 14, 52, 13, 52, 51, 14, 15, 53, 14, 53, 52, 15, 16, 54, 15, 54, 53, 16, 17, 55, 16, 55, 54, 19, 20, 58, 19, 58, 57, 20, 21, 59, 20, 59, 58, 21, 22, 60, 21, 60, 59, 22, 23, 61, 22, 61, 60, 25, 26, 64, 25, 64, 63, 26, 27, 65, 26, 65, 64, 27, 28, 66, 27, 66, 65, 28, 29, 67, 28, 67, 66, 29, 30, 68, 29, 68, 67, 30, 31, 69, 30, 69, 68, 31, 32, 70, 31, 70, 69, 32, 33, 71, 32, 71, 70, 33, 34, 72, 33, 72, 71, 34, 35, 73, 34, 73, 72, 35, 36, 74, 35, 74, 73, 37, 0, 36, 76, 36, 0, 76, 74, 36, 76, 38, 74, 75, 74, 38, 4, 5, 6, 4, 6, 77, 77, 6, 44, 42, 77, 44, 42, 44, 43, 18, 19, 17, 17, 19, 78, 17, 78, 55, 55, 78, 57, 56, 55, 57, 25, 23, 24, 25, 79, 23, 79, 25, 63, 61, 79, 63, 61, 63, 62], [ 0, 2, 1, 1, 2, 38, 2, 39, 38, 0, 3, 2, 2, 3, 39, 3, 40, 39, 0, 4, 3, 3, 4, 40, 4, 41, 40, 0, 5, 4, 4, 5, 41, 5, 42, 41, 0, 6, 5, 5, 6, 42, 6, 43, 42, 0, 7, 6, 6, 7, 43, 7, 44, 43, 0, 8, 7, 7, 8, 44, 8, 45, 44, 0, 9, 8, 8, 9, 45, 9, 46, 45, 0, 10, 9, 9, 10, 46, 10, 47, 46, 0, 11, 10, 10, 11, 47, 11, 48, 47, 0, 12, 11, 11, 12, 48, 12, 49, 48, 0, 13, 12, 12, 13, 49, 13, 50, 49, 0, 14, 13, 13, 14, 50, 14, 51, 50, 0, 15, 14, 14, 15, 51, 15, 52, 51, 0, 16, 15, 15, 16, 52, 16, 53, 52, 0, 17, 16, 16, 17, 53, 17, 54, 53, 0, 18, 17, 17, 18, 54, 18, 55, 54, 0, 19, 18, 18, 19, 55, 19, 56, 55, 0, 20, 19, 19, 20, 56, 20, 57, 56, 0, 21, 20, 20, 21, 57, 21, 58, 57, 0, 22, 21, 21, 22, 58, 22, 59, 58, 0, 23, 22, 22, 23, 59, 23, 60, 59, 0, 24, 23, 23, 24, 60, 24, 61, 60, 0, 25, 24, 24, 25, 61, 25, 62, 61, 0, 26, 25, 25, 26, 62, 26, 63, 62, 0, 27, 26, 26, 27, 63, 27, 64, 63, 0, 28, 27, 27, 28, 64, 28, 65, 64, 0, 29, 28, 28, 29, 65, 29, 66, 65, 0, 30, 29, 29, 30, 66, 30, 67, 66, 0, 31, 30, 30, 31, 67, 31, 68, 67, 0, 32, 31, 31, 32, 68, 32, 69, 68, 0, 33, 32, 32, 33, 69, 33, 70, 69, 0, 34, 33, 33, 34, 70, 34, 71, 70, 0, 35, 34, 34, 35, 71, 35, 72, 71, 0, 36, 35, 35, 36, 72, 36, 73, 72, 0, 1, 36, 36, 1, 73, 1, 38, 73, 37, 38, 39, 37, 39, 40, 37, 40, 41, 37, 41, 42, 37, 42, 43, 37, 43, 44, 37, 44, 45, 37, 45, 46, 37, 46, 47, 37, 47, 48, 37, 48, 49, 37, 49, 50, 37, 50, 51, 37, 51, 52, 37, 52, 53, 37, 53, 54, 37, 54, 55, 37, 55, 56, 37, 56, 57, 37, 57, 58, 37, 58, 59, 37, 59, 60, 37, 60, 61, 37, 61, 62, 37, 62, 63, 37, 63, 64, 37, 64, 65, 37, 65, 66, 37, 66, 67, 37, 67, 68, 37, 68, 69, 37, 69, 70, 37, 70, 71, 37, 71, 72, 37, 72, 73, 37, 73, 38], [ 0, 2, 1, 1, 2, 38, 2, 39, 38, 0, 3, 2, 2, 3, 39, 3, 40, 39, 0, 4, 3, 3, 4, 40, 4, 41, 40, 0, 5, 4, 4, 5, 41, 5, 42, 41, 0, 6, 5, 5, 6, 42, 6, 43, 42, 0, 7, 6, 6, 7, 43, 7, 44, 43, 0, 8, 7, 7, 8, 44, 8, 45, 44, 0, 9, 8, 8, 9, 45, 9, 46, 45, 0, 10, 9, 9, 10, 46, 10, 47, 46, 0, 11, 10, 10, 11, 47, 11, 48, 47, 0, 12, 11, 11, 12, 48, 12, 49, 48, 0, 13, 12, 12, 13, 49, 13, 50, 49, 0, 14, 13, 13, 14, 50, 14, 51, 50, 0, 15, 14, 14, 15, 51, 15, 52, 51, 0, 16, 15, 15, 16, 52, 16, 53, 52, 0, 17, 16, 16, 17, 53, 17, 54, 53, 0, 18, 17, 17, 18, 54, 18, 55, 54, 0, 19, 18, 18, 19, 55, 19, 56, 55, 0, 20, 19, 19, 20, 56, 20, 57, 56, 0, 21, 20, 20, 21, 57, 21, 58, 57, 0, 22, 21, 21, 22, 58, 22, 59, 58, 0, 23, 22, 22, 23, 59, 23, 60, 59, 0, 24, 23, 23, 24, 60, 24, 61, 60, 0, 25, 24, 24, 25, 61, 25, 62, 61, 0, 26, 25, 25, 26, 62, 26, 63, 62, 0, 27, 26, 26, 27, 63, 27, 64, 63, 0, 28, 27, 27, 28, 64, 28, 65, 64, 0, 29, 28, 28, 29, 65, 29, 66, 65, 0, 30, 29, 29, 30, 66, 30, 67, 66, 0, 31, 30, 30, 31, 67, 31, 68, 67, 0, 32, 31, 31, 32, 68, 32, 69, 68, 0, 33, 32, 32, 33, 69, 33, 70, 69, 0, 34, 33, 33, 34, 70, 34, 71, 70, 0, 35, 34, 34, 35, 71, 35, 72, 71, 0, 36, 35, 35, 36, 72, 36, 73, 72, 0, 1, 36, 36, 1, 73, 1, 38, 73, 37, 38, 39, 37, 39, 40, 37, 40, 41, 37, 41, 42, 37, 42, 43, 37, 43, 44, 37, 44, 45, 37, 45, 46, 37, 46, 47, 37, 47, 48, 37, 48, 49, 37, 49, 50, 37, 50, 51, 37, 51, 52, 37, 52, 53, 37, 53, 54, 37, 54, 55, 37, 55, 56, 37, 56, 57, 37, 57, 58, 37, 58, 59, 37, 59, 60, 37, 60, 61, 37, 61, 62, 37, 62, 63, 37, 63, 64, 37, 64, 65, 37, 65, 66, 37, 66, 67, 37, 67, 68, 37, 68, 69, 37, 69, 70, 37, 70, 71, 37, 71, 72, 37, 72, 73, 37, 73, 38], [ 0, 2, 1, 1, 2, 38, 2, 39, 38, 0, 3, 2, 2, 3, 39, 3, 40, 39, 0, 4, 3, 3, 4, 40, 4, 41, 40, 0, 5, 4, 4, 5, 41, 5, 42, 41, 0, 6, 5, 5, 6, 42, 6, 43, 42, 0, 7, 6, 6, 7, 43, 7, 44, 43, 0, 8, 7, 7, 8, 44, 8, 45, 44, 0, 9, 8, 8, 9, 45, 9, 46, 45, 0, 10, 9, 9, 10, 46, 10, 47, 46, 0, 11, 10, 10, 11, 47, 11, 48, 47, 0, 12, 11, 11, 12, 48, 12, 49, 48, 0, 13, 12, 12, 13, 49, 13, 50, 49, 0, 14, 13, 13, 14, 50, 14, 51, 50, 0, 15, 14, 14, 15, 51, 15, 52, 51, 0, 16, 15, 15, 16, 52, 16, 53, 52, 0, 17, 16, 16, 17, 53, 17, 54, 53, 0, 18, 17, 17, 18, 54, 18, 55, 54, 0, 19, 18, 18, 19, 55, 19, 56, 55, 0, 20, 19, 19, 20, 56, 20, 57, 56, 0, 21, 20, 20, 21, 57, 21, 58, 57, 0, 22, 21, 21, 22, 58, 22, 59, 58, 0, 23, 22, 22, 23, 59, 23, 60, 59, 0, 24, 23, 23, 24, 60, 24, 61, 60, 0, 25, 24, 24, 25, 61, 25, 62, 61, 0, 26, 25, 25, 26, 62, 26, 63, 62, 0, 27, 26, 26, 27, 63, 27, 64, 63, 0, 28, 27, 27, 28, 64, 28, 65, 64, 0, 29, 28, 28, 29, 65, 29, 66, 65, 0, 30, 29, 29, 30, 66, 30, 67, 66, 0, 31, 30, 30, 31, 67, 31, 68, 67, 0, 32, 31, 31, 32, 68, 32, 69, 68, 0, 33, 32, 32, 33, 69, 33, 70, 69, 0, 34, 33, 33, 34, 70, 34, 71, 70, 0, 35, 34, 34, 35, 71, 35, 72, 71, 0, 36, 35, 35, 36, 72, 36, 73, 72, 0, 1, 36, 36, 1, 73, 1, 38, 73, 37, 38, 39, 37, 39, 40, 37, 40, 41, 37, 41, 42, 37, 42, 43, 37, 43, 44, 37, 44, 45, 37, 45, 46, 37, 46, 47, 37, 47, 48, 37, 48, 49, 37, 49, 50, 37, 50, 51, 37, 51, 52, 37, 52, 53, 37, 53, 54, 37, 54, 55, 37, 55, 56, 37, 56, 57, 37, 57, 58, 37, 58, 59, 37, 59, 60, 37, 60, 61, 37, 61, 62, 37, 62, 63, 37, 63, 64, 37, 64, 65, 37, 65, 66, 37, 66, 67, 37, 67, 68, 37, 68, 69, 37, 69, 70, 37, 70, 71, 37, 71, 72, 37, 72, 73, 37, 73, 38], [ 0, 2, 1, 1, 2, 38, 2, 39, 38, 0, 3, 2, 2, 3, 39, 3, 40, 39, 0, 4, 3, 3, 4, 40, 4, 41, 40, 0, 5, 4, 4, 5, 41, 5, 42, 41, 0, 6, 5, 5, 6, 42, 6, 43, 42, 0, 7, 6, 6, 7, 43, 7, 44, 43, 0, 8, 7, 7, 8, 44, 8, 45, 44, 0, 9, 8, 8, 9, 45, 9, 46, 45, 0, 10, 9, 9, 10, 46, 10, 47, 46, 0, 11, 10, 10, 11, 47, 11, 48, 47, 0, 12, 11, 11, 12, 48, 12, 49, 48, 0, 13, 12, 12, 13, 49, 13, 50, 49, 0, 14, 13, 13, 14, 50, 14, 51, 50, 0, 15, 14, 14, 15, 51, 15, 52, 51, 0, 16, 15, 15, 16, 52, 16, 53, 52, 0, 17, 16, 16, 17, 53, 17, 54, 53, 0, 18, 17, 17, 18, 54, 18, 55, 54, 0, 19, 18, 18, 19, 55, 19, 56, 55, 0, 20, 19, 19, 20, 56, 20, 57, 56, 0, 21, 20, 20, 21, 57, 21, 58, 57, 0, 22, 21, 21, 22, 58, 22, 59, 58, 0, 23, 22, 22, 23, 59, 23, 60, 59, 0, 24, 23, 23, 24, 60, 24, 61, 60, 0, 25, 24, 24, 25, 61, 25, 62, 61, 0, 26, 25, 25, 26, 62, 26, 63, 62, 0, 27, 26, 26, 27, 63, 27, 64, 63, 0, 28, 27, 27, 28, 64, 28, 65, 64, 0, 29, 28, 28, 29, 65, 29, 66, 65, 0, 30, 29, 29, 30, 66, 30, 67, 66, 0, 31, 30, 30, 31, 67, 31, 68, 67, 0, 32, 31, 31, 32, 68, 32, 69, 68, 0, 33, 32, 32, 33, 69, 33, 70, 69, 0, 34, 33, 33, 34, 70, 34, 71, 70, 0, 35, 34, 34, 35, 71, 35, 72, 71, 0, 36, 35, 35, 36, 72, 36, 73, 72, 0, 1, 36, 36, 1, 73, 1, 38, 73, 37, 38, 39, 37, 39, 40, 37, 40, 41, 37, 41, 42, 37, 42, 43, 37, 43, 44, 37, 44, 45, 37, 45, 46, 37, 46, 47, 37, 47, 48, 37, 48, 49, 37, 49, 50, 37, 50, 51, 37, 51, 52, 37, 52, 53, 37, 53, 54, 37, 54, 55, 37, 55, 56, 37, 56, 57, 37, 57, 58, 37, 58, 59, 37, 59, 60, 37, 60, 61, 37, 61, 62, 37, 62, 63, 37, 63, 64, 37, 64, 65, 37, 65, 66, 37, 66, 67, 37, 67, 68, 37, 68, 69, 37, 69, 70, 37, 70, 71, 37, 71, 72, 37, 72, 73, 37, 73, 38], [ 0, 2, 1, 1, 2, 38, 2, 39, 38, 0, 3, 2, 2, 3, 39, 3, 40, 39, 0, 4, 3, 3, 4, 40, 4, 41, 40, 0, 5, 4, 4, 5, 41, 5, 42, 41, 0, 6, 5, 5, 6, 42, 6, 43, 42, 0, 7, 6, 6, 7, 43, 7, 44, 43, 0, 8, 7, 7, 8, 44, 8, 45, 44, 0, 9, 8, 8, 9, 45, 9, 46, 45, 0, 10, 9, 9, 10, 46, 10, 47, 46, 0, 11, 10, 10, 11, 47, 11, 48, 47, 0, 12, 11, 11, 12, 48, 12, 49, 48, 0, 13, 12, 12, 13, 49, 13, 50, 49, 0, 14, 13, 13, 14, 50, 14, 51, 50, 0, 15, 14, 14, 15, 51, 15, 52, 51, 0, 16, 15, 15, 16, 52, 16, 53, 52, 0, 17, 16, 16, 17, 53, 17, 54, 53, 0, 18, 17, 17, 18, 54, 18, 55, 54, 0, 19, 18, 18, 19, 55, 19, 56, 55, 0, 20, 19, 19, 20, 56, 20, 57, 56, 0, 21, 20, 20, 21, 57, 21, 58, 57, 0, 22, 21, 21, 22, 58, 22, 59, 58, 0, 23, 22, 22, 23, 59, 23, 60, 59, 0, 24, 23, 23, 24, 60, 24, 61, 60, 0, 25, 24, 24, 25, 61, 25, 62, 61, 0, 26, 25, 25, 26, 62, 26, 63, 62, 0, 27, 26, 26, 27, 63, 27, 64, 63, 0, 28, 27, 27, 28, 64, 28, 65, 64, 0, 29, 28, 28, 29, 65, 29, 66, 65, 0, 30, 29, 29, 30, 66, 30, 67, 66, 0, 31, 30, 30, 31, 67, 31, 68, 67, 0, 32, 31, 31, 32, 68, 32, 69, 68, 0, 33, 32, 32, 33, 69, 33, 70, 69, 0, 34, 33, 33, 34, 70, 34, 71, 70, 0, 35, 34, 34, 35, 71, 35, 72, 71, 0, 36, 35, 35, 36, 72, 36, 73, 72, 0, 1, 36, 36, 1, 73, 1, 38, 73, 37, 38, 39, 37, 39, 40, 37, 40, 41, 37, 41, 42, 37, 42, 43, 37, 43, 44, 37, 44, 45, 37, 45, 46, 37, 46, 47, 37, 47, 48, 37, 48, 49, 37, 49, 50, 37, 50, 51, 37, 51, 52, 37, 52, 53, 37, 53, 54, 37, 54, 55, 37, 55, 56, 37, 56, 57, 37, 57, 58, 37, 58, 59, 37, 59, 60, 37, 60, 61, 37, 61, 62, 37, 62, 63, 37, 63, 64, 37, 64, 65, 37, 65, 66, 37, 66, 67, 37, 67, 68, 37, 68, 69, 37, 69, 70, 37, 70, 71, 37, 71, 72, 37, 72, 73, 37, 73, 38], [ 2, 1, 0, 2, 4, 3, 2, 5, 4, 2, 0, 49, 2, 49, 5, 49, 48, 5, 5, 48, 6, 27, 26, 25, 27, 29, 28, 27, 30, 29, 27, 25, 24, 27, 24, 30, 24, 23, 30, 30, 23, 31, 52, 50, 51, 52, 53, 54, 52, 54, 55, 52, 99, 50, 52, 55, 99, 99, 55, 98, 55, 56, 98, 77, 75, 76, 77, 78, 79, 77, 79, 80, 77, 74, 75, 77, 80, 74, 74, 80, 73, 80, 81, 73, 6, 48, 7, 48, 47, 7, 56, 57, 98, 98, 57, 97, 7, 47, 8, 47, 46, 8, 57, 58, 97, 97, 58, 96, 8, 46, 9, 46, 45, 9, 58, 59, 96, 96, 59, 95, 9, 45, 10, 45, 44, 10, 59, 60, 95, 95, 60, 94, 10, 44, 11, 44, 43, 11, 60, 61, 94, 94, 61, 93, 11, 43, 12, 43, 42, 12, 61, 62, 93, 93, 62, 92, 12, 42, 13, 42, 41, 13, 62, 63, 92, 92, 63, 91, 13, 41, 14, 41, 40, 14, 63, 64, 91, 91, 64, 90, 14, 40, 15, 40, 39, 15, 64, 65, 90, 90, 65, 89, 15, 39, 16, 39, 38, 16, 65, 66, 89, 89, 66, 88, 16, 38, 17, 38, 37, 17, 66, 67, 88, 88, 67, 87, 17, 37, 18, 37, 36, 18, 67, 68, 87, 87, 68, 86, 18, 36, 19, 36, 35, 19, 68, 69, 86, 86, 69, 85, 19, 35, 20, 35, 34, 20, 69, 70, 85, 85, 70, 84, 20, 34, 21, 34, 33, 21, 70, 71, 84, 84, 71, 83, 21, 33, 22, 33, 32, 22, 71, 72, 83, 83, 72, 82, 22, 32, 23, 32, 31, 23, 72, 73, 82, 82, 73, 81, 0, 1, 51, 0, 51, 50, 1, 2, 52, 1, 52, 51, 2, 3, 53, 2, 53, 52, 3, 4, 54, 3, 54, 53, 6, 7, 57, 6, 57, 56, 7, 8, 58, 7, 58, 57, 8, 9, 59, 8, 59, 58, 9, 10, 60, 9, 60, 59, 10, 11, 61, 10, 61, 60, 11, 12, 62, 11, 62, 61, 12, 13, 63, 12, 63, 62, 13, 14, 64, 13, 64, 63, 14, 15, 65, 14, 65, 64, 15, 16, 66, 15, 66, 65, 16, 17, 67, 16, 67, 66, 17, 18, 68, 17, 68, 67, 18, 19, 69, 18, 69, 68, 19, 20, 70, 19, 70, 69, 20, 21, 71, 20, 71, 70, 21, 22, 72, 21, 72, 71, 22, 23, 73, 22, 73, 72, 25, 26, 76, 25, 76, 75, 26, 27, 77, 26, 77, 76, 27, 28, 78, 27, 78, 77, 28, 29, 79, 28, 79, 78, 31, 32, 82, 31, 82, 81, 32, 33, 83, 32, 83, 82, 33, 34, 84, 33, 84, 83, 34, 35, 85, 34, 85, 84, 35, 36, 86, 35, 86, 85, 36, 37, 87, 36, 87, 86, 37, 38, 88, 37, 88, 87, 38, 39, 89, 38, 89, 88, 39, 40, 90, 39, 90, 89, 40, 41, 91, 40, 91, 90, 41, 42, 92, 41, 92, 91, 42, 43, 93, 42, 93, 92, 43, 44, 94, 43, 94, 93, 44, 45, 95, 44, 95, 94, 45, 46, 96, 45, 96, 95, 46, 47, 97, 46, 97, 96, 47, 48, 98, 47, 98, 97, 49, 0, 48, 100, 48, 0, 100, 98, 48, 100, 50, 98, 99, 98, 50, 4, 5, 6, 4, 6, 101, 101, 6, 56, 54, 101, 56, 54, 56, 55, 24, 25, 23, 23, 25, 102, 23, 102, 73, 73, 102, 75, 74, 73, 75, 31, 29, 30, 31, 103, 29, 103, 31, 81, 79, 103, 81, 79, 81, 80], [ 0, 2, 1, 1, 2, 38, 2, 39, 38, 0, 3, 2, 2, 3, 39, 3, 40, 39, 0, 4, 3, 3, 4, 40, 4, 41, 40, 0, 5, 4, 4, 5, 41, 5, 42, 41, 0, 6, 5, 5, 6, 42, 6, 43, 42, 0, 7, 6, 6, 7, 43, 7, 44, 43, 0, 8, 7, 7, 8, 44, 8, 45, 44, 0, 9, 8, 8, 9, 45, 9, 46, 45, 0, 10, 9, 9, 10, 46, 10, 47, 46, 0, 11, 10, 10, 11, 47, 11, 48, 47, 0, 12, 11, 11, 12, 48, 12, 49, 48, 0, 13, 12, 12, 13, 49, 13, 50, 49, 0, 14, 13, 13, 14, 50, 14, 51, 50, 0, 15, 14, 14, 15, 51, 15, 52, 51, 0, 16, 15, 15, 16, 52, 16, 53, 52, 0, 17, 16, 16, 17, 53, 17, 54, 53, 0, 18, 17, 17, 18, 54, 18, 55, 54, 0, 19, 18, 18, 19, 55, 19, 56, 55, 0, 20, 19, 19, 20, 56, 20, 57, 56, 0, 21, 20, 20, 21, 57, 21, 58, 57, 0, 22, 21, 21, 22, 58, 22, 59, 58, 0, 23, 22, 22, 23, 59, 23, 60, 59, 0, 24, 23, 23, 24, 60, 24, 61, 60, 0, 25, 24, 24, 25, 61, 25, 62, 61, 0, 26, 25, 25, 26, 62, 26, 63, 62, 0, 27, 26, 26, 27, 63, 27, 64, 63, 0, 28, 27, 27, 28, 64, 28, 65, 64, 0, 29, 28, 28, 29, 65, 29, 66, 65, 0, 30, 29, 29, 30, 66, 30, 67, 66, 0, 31, 30, 30, 31, 67, 31, 68, 67, 0, 32, 31, 31, 32, 68, 32, 69, 68, 0, 33, 32, 32, 33, 69, 33, 70, 69, 0, 34, 33, 33, 34, 70, 34, 71, 70, 0, 35, 34, 34, 35, 71, 35, 72, 71, 0, 36, 35, 35, 36, 72, 36, 73, 72, 0, 1, 36, 36, 1, 73, 1, 38, 73, 37, 38, 39, 37, 39, 40, 37, 40, 41, 37, 41, 42, 37, 42, 43, 37, 43, 44, 37, 44, 45, 37, 45, 46, 37, 46, 47, 37, 47, 48, 37, 48, 49, 37, 49, 50, 37, 50, 51, 37, 51, 52, 37, 52, 53, 37, 53, 54, 37, 54, 55, 37, 55, 56, 37, 56, 57, 37, 57, 58, 37, 58, 59, 37, 59, 60, 37, 60, 61, 37, 61, 62, 37, 62, 63, 37, 63, 64, 37, 64, 65, 37, 65, 66, 37, 66, 67, 37, 67, 68, 37, 68, 69, 37, 69, 70, 37, 70, 71, 37, 71, 72, 37, 72, 73, 37, 73, 38], [ 0, 2, 1, 1, 2, 38, 2, 39, 38, 0, 3, 2, 2, 3, 39, 3, 40, 39, 0, 4, 3, 3, 4, 40, 4, 41, 40, 0, 5, 4, 4, 5, 41, 5, 42, 41, 0, 6, 5, 5, 6, 42, 6, 43, 42, 0, 7, 6, 6, 7, 43, 7, 44, 43, 0, 8, 7, 7, 8, 44, 8, 45, 44, 0, 9, 8, 8, 9, 45, 9, 46, 45, 0, 10, 9, 9, 10, 46, 10, 47, 46, 0, 11, 10, 10, 11, 47, 11, 48, 47, 0, 12, 11, 11, 12, 48, 12, 49, 48, 0, 13, 12, 12, 13, 49, 13, 50, 49, 0, 14, 13, 13, 14, 50, 14, 51, 50, 0, 15, 14, 14, 15, 51, 15, 52, 51, 0, 16, 15, 15, 16, 52, 16, 53, 52, 0, 17, 16, 16, 17, 53, 17, 54, 53, 0, 18, 17, 17, 18, 54, 18, 55, 54, 0, 19, 18, 18, 19, 55, 19, 56, 55, 0, 20, 19, 19, 20, 56, 20, 57, 56, 0, 21, 20, 20, 21, 57, 21, 58, 57, 0, 22, 21, 21, 22, 58, 22, 59, 58, 0, 23, 22, 22, 23, 59, 23, 60, 59, 0, 24, 23, 23, 24, 60, 24, 61, 60, 0, 25, 24, 24, 25, 61, 25, 62, 61, 0, 26, 25, 25, 26, 62, 26, 63, 62, 0, 27, 26, 26, 27, 63, 27, 64, 63, 0, 28, 27, 27, 28, 64, 28, 65, 64, 0, 29, 28, 28, 29, 65, 29, 66, 65, 0, 30, 29, 29, 30, 66, 30, 67, 66, 0, 31, 30, 30, 31, 67, 31, 68, 67, 0, 32, 31, 31, 32, 68, 32, 69, 68, 0, 33, 32, 32, 33, 69, 33, 70, 69, 0, 34, 33, 33, 34, 70, 34, 71, 70, 0, 35, 34, 34, 35, 71, 35, 72, 71, 0, 36, 35, 35, 36, 72, 36, 73, 72, 0, 1, 36, 36, 1, 73, 1, 38, 73, 37, 38, 39, 37, 39, 40, 37, 40, 41, 37, 41, 42, 37, 42, 43, 37, 43, 44, 37, 44, 45, 37, 45, 46, 37, 46, 47, 37, 47, 48, 37, 48, 49, 37, 49, 50, 37, 50, 51, 37, 51, 52, 37, 52, 53, 37, 53, 54, 37, 54, 55, 37, 55, 56, 37, 56, 57, 37, 57, 58, 37, 58, 59, 37, 59, 60, 37, 60, 61, 37, 61, 62, 37, 62, 63, 37, 63, 64, 37, 64, 65, 37, 65, 66, 37, 66, 67, 37, 67, 68, 37, 68, 69, 37, 69, 70, 37, 70, 71, 37, 71, 72, 37, 72, 73, 37, 73, 38], [ 0, 2, 1, 1, 2, 38, 2, 39, 38, 0, 3, 2, 2, 3, 39, 3, 40, 39, 0, 4, 3, 3, 4, 40, 4, 41, 40, 0, 5, 4, 4, 5, 41, 5, 42, 41, 0, 6, 5, 5, 6, 42, 6, 43, 42, 0, 7, 6, 6, 7, 43, 7, 44, 43, 0, 8, 7, 7, 8, 44, 8, 45, 44, 0, 9, 8, 8, 9, 45, 9, 46, 45, 0, 10, 9, 9, 10, 46, 10, 47, 46, 0, 11, 10, 10, 11, 47, 11, 48, 47, 0, 12, 11, 11, 12, 48, 12, 49, 48, 0, 13, 12, 12, 13, 49, 13, 50, 49, 0, 14, 13, 13, 14, 50, 14, 51, 50, 0, 15, 14, 14, 15, 51, 15, 52, 51, 0, 16, 15, 15, 16, 52, 16, 53, 52, 0, 17, 16, 16, 17, 53, 17, 54, 53, 0, 18, 17, 17, 18, 54, 18, 55, 54, 0, 19, 18, 18, 19, 55, 19, 56, 55, 0, 20, 19, 19, 20, 56, 20, 57, 56, 0, 21, 20, 20, 21, 57, 21, 58, 57, 0, 22, 21, 21, 22, 58, 22, 59, 58, 0, 23, 22, 22, 23, 59, 23, 60, 59, 0, 24, 23, 23, 24, 60, 24, 61, 60, 0, 25, 24, 24, 25, 61, 25, 62, 61, 0, 26, 25, 25, 26, 62, 26, 63, 62, 0, 27, 26, 26, 27, 63, 27, 64, 63, 0, 28, 27, 27, 28, 64, 28, 65, 64, 0, 29, 28, 28, 29, 65, 29, 66, 65, 0, 30, 29, 29, 30, 66, 30, 67, 66, 0, 31, 30, 30, 31, 67, 31, 68, 67, 0, 32, 31, 31, 32, 68, 32, 69, 68, 0, 33, 32, 32, 33, 69, 33, 70, 69, 0, 34, 33, 33, 34, 70, 34, 71, 70, 0, 35, 34, 34, 35, 71, 35, 72, 71, 0, 36, 35, 35, 36, 72, 36, 73, 72, 0, 1, 36, 36, 1, 73, 1, 38, 73, 37, 38, 39, 37, 39, 40, 37, 40, 41, 37, 41, 42, 37, 42, 43, 37, 43, 44, 37, 44, 45, 37, 45, 46, 37, 46, 47, 37, 47, 48, 37, 48, 49, 37, 49, 50, 37, 50, 51, 37, 51, 52, 37, 52, 53, 37, 53, 54, 37, 54, 55, 37, 55, 56, 37, 56, 57, 37, 57, 58, 37, 58, 59, 37, 59, 60, 37, 60, 61, 37, 61, 62, 37, 62, 63, 37, 63, 64, 37, 64, 65, 37, 65, 66, 37, 66, 67, 37, 67, 68, 37, 68, 69, 37, 69, 70, 37, 70, 71, 37, 71, 72, 37, 72, 73, 37, 73, 38]]
decoration = [[ 0.000000, 0.000000, 32.632000, 0.000000, 0.000000, 33.000000], [ 0.000000, 0.000000, 33.202000, 0.000000, 0.000000, 33.570000], [ 0.000000, 0.000000, 33.772000, 0.000000, 0.000000, 34.140000], [ 0.000000, 0.000000, 34.452500, 0.000000, 0.000000, 34.462500], [ 0.000000, 0.000000, 34.452500, 0.000000, 0.000000, 34.462500], [ 0.000000, 0.000000, 36.192100, 0.000000, 0.000000, 36.202100], [ 0.000000, 0.000000, 36.192100, 0.000000, 0.000000, 36.202100], [ 0.000000, 0.000000, 36.287600, 0.000000, 0.000000, 36.297600], [ 0.000000, 0.000000, 36.287600, 0.000000, 0.000000, 36.297600], [ 0.000000, 0.000000, 37.366700, 0.000000, 0.000000, 37.376700], [ 0.000000, 0.000000, 37.366700, 0.000000, 0.000000, 37.376700], [ 0.000000, 0.000000, 37.650200, 0.000000, 0.000000, 37.970200], [ 0.000000, 0.000000, 38.280200, 0.000000, 0.000000, 38.600200], [ -0.393380, 0.000000, 39.036200, 0.393380, 0.000000, 39.036200, -0.924534, 0.000000, 40.021871, -0.491878, 0.000000, 40.678985, 0.000000, 0.000000, 39.036200, 0.000000, 0.000000, 39.036200, -0.708206, 0.000000, 40.350428, -0.708206, 0.000000, 40.350428], [ -1.036641, 0.000000, 40.566676, -1.044993, 0.000000, 40.572175], [ -1.036641, 0.000000, 40.566676, -1.044993, 0.000000, 40.572175], [ -1.263569, 0.000000, 40.716089, -1.530839, 0.000000, 40.892064], [ -1.814812, 0.000000, 41.079038, -2.082081, 0.000000, 41.255013], [ -2.653454, 0.000000, 41.160054, -2.220642, 0.000000, 41.817405, -3.768777, 0.000000, 41.258684, -4.079513, 0.000000, 41.981788, -2.437048, 0.000000, 41.488730, -2.437048, 0.000000, 41.488730, -3.924145, 0.000000, 41.620236, -3.924145, 0.000000, 41.620236], [ -4.285422, 0.000000, 41.464986, -4.294610, 0.000000, 41.461038], [ -4.285422, 0.000000, 41.464986, -4.294610, 0.000000, 41.461038], [ -4.723212, 0.000000, 41.276857, -5.017215, 0.000000, 41.150516], [ -5.145842, 0.000000, 41.095242, -5.439845, 0.000000, 40.968901], [ -5.568472, 0.000000, 40.913627, -5.862475, 0.000000, 40.787286], [ -7.423303, 0.000000, 40.540911, -7.115444, -0.000000, 39.824502, -8.947706, 0.000000, 41.126758, -9.656116, 0.000000, 40.800915, -7.269373, -0.000000, 40.182706, -7.269373, -0.000000, 40.182706, -9.301911, 0.000000, 40.963837, -9.301912, 0.000000, 40.963837], [ -10.729840, 0.000000, 44.068266, -10.734019, 0.000000, 44.077351], [ -10.729840, 0.000000, 44.068266, -10.734019, 0.000000, 44.077351], [ -11.363260, 0.000000, 45.445371, -11.367439, 0.000000, 45.454456]]

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
        if math.abs(dot(u, ref) - 1.0) < 1e-12:
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
