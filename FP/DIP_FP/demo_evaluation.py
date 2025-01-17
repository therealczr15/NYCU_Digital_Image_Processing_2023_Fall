import cv2
import numpy as np

def calculate_iou(mask1, mask2):
    intersection = np.logical_and(mask1, mask2)
    union = np.logical_or(mask1, mask2)
    iou = np.sum(intersection) / np.sum(union)
    return iou

# 假設你有10對mask的檔案名稱，這裡只是一個例子，請替換成你的實際檔案名稱
mask_filenames_1 = ["predict_mask/ADE_train_00000591.png", "predict_mask/ADE_train_00000650.png", "predict_mask/ADE_train_00000863.png", "predict_mask/ADE_train_00000941.png", "predict_mask/ADE_train_00000964.png"]

mask_filenames_2 = ["test_dataset/mask/ADE_train_00000591.png", "test_dataset/mask/ADE_train_00000650.png", "test_dataset/mask/ADE_train_00000863.png", "test_dataset/mask/ADE_train_00000941.png", "test_dataset/mask/ADE_train_00000964.png"]

total_iou = 0

for i in range(5):  # 修改為實際的mask數量
    # 讀取兩個mask
    mask1 = cv2.imread(mask_filenames_1[i], cv2.IMREAD_GRAYSCALE)
    mask2 = cv2.imread(mask_filenames_2[i], cv2.IMREAD_GRAYSCALE)

    # 將mask二值化 (假設是二值化的黑白圖片)
    _, mask1 = cv2.threshold(mask1, 128, 255, cv2.THRESH_BINARY)
    _, mask2 = cv2.threshold(mask2, 128, 255, cv2.THRESH_BINARY)

    # 計算IOU並加總
    iou = calculate_iou(mask1, mask2)
    total_iou += iou
    print(f"IoU for pair {i + 1}: {iou}")

# 計算平均IoU
average_iou = total_iou / len(mask_filenames_1)
print(f"Average IoU: {average_iou}")
