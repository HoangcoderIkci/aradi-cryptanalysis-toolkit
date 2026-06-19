# Tấn công cube lên ARADI — giải thích dễ hiểu

> **Phạm vi (đọc trước):** Bài viết mô tả việc **kiểm chứng độc lập** một cube-distinguisher
> của ARADI và một **tấn công thực nghiệm lên biến thể 6 vòng ĐÃ SỬA key schedule**. Nó
> **không phá** ARADI đầy đủ (16 vòng vẫn ~2¹⁴¹, thực tế bất khả thi) và **không** khôi phục
> master key 256-bit (chỉ phục hồi cặp 32-bit `(u, ℓ)`). Xem [`../SCOPE.md`](../SCOPE.md).

Mã nguồn tái lập: [repo này](../README.md) — `make test-c`, `make reproduce`, `make milp`.

---

## 1. ARADI là gì

ARADI là block cipher do **NSA công bố 2024** (ePrint 2024/1240), thiết kế cho **mã hoá bộ nhớ
(memory encryption)** nên ưu tiên **độ trễ thấp**:

- Khối **128 bit** (4 từ 32-bit: W, X, Y, Z), khoá **256 bit**, **16 vòng**.
- Mỗi vòng = XOR khoá vòng → **S-box** → **lớp tuyến tính** (linear layer).
- S-box dựng từ **4 cổng Toffoli**; cả khối S-box có **bậc đại số = 3** — *thấp*. Bậc thấp chính
  là thứ khiến ARADI thành mục tiêu tự nhiên cho **kubическая (cube) attack**.

## 2. Cube attack — trực giác, không cần đại số nặng

Mọi bit đầu ra của cipher là một đa thức Boolean theo bit đầu vào. Ý tưởng:

1. Chọn **13 bit** của bản rõ làm "cube", giữ các bit còn lại cố định.
2. Cho 13 bit đó chạy hết **2¹³ = 8192** tổ hợp, **cộng XOR** tất cả bản mã thu được.
3. Phép cộng trên một cube chiều `k` **triệt tiêu mọi đơn thức bậc < k**. Vì bậc của ARADI thấp,
   sau một số vòng tổng này lộ ra **cấu trúc** thay vì ngẫu nhiên — đó là **distinguisher**.

Điểm mấu chốt: ta không phân tích **một** bản mã, mà phân tích **tổng theo một tập** bản mã.

## 3. AABB distinguisher (sau 5 vòng)

Sau **5 vòng**, tổng-cube ở các từ trạng thái **X** và **Z** xuất hiện cấu trúc byte **trong cùng
một từ (intra-word)**: hai byte đầu bằng nhau, hai byte cuối bằng nhau — dạng **A-A-B-B**.

- Với cipher ngẫu nhiên, xác suất trùng kiểu này cỡ **2⁻³²**.
- Với ARADI, nó xảy ra **tất định (deterministic)**.

Mình **kiểm chứng độc lập** tính chất này bằng **MILP** (mô hình lan truyền bậc, giải bằng
PuLP + CBC — solver miễn phí), và bằng cài đặt C **khớp test-vector chính thức của NSA**. Nghĩa
là AABB là **tính chất cấu trúc của cipher**, không phải lỗi cài đặt. Chi tiết:
[`../python/aradi_milp.py`](../python/aradi_milp.py), [`../python/run_classification.py`](../python/run_classification.py).

## 4. ARADI đầy đủ thì sao? — Vì sao KHÔNG phá được

Từ distinguisher có thể dựng một tấn công khôi phục khoá **2 giai đoạn**:

- **Giai đoạn 1:** dò khoá vòng cuối `RK₅` (128 bit) + lọc bằng AABB.
- **Giai đoạn 2:** dò `RK₄` (chỉ 64 bit hiệu dụng nhờ quan sát của Bellini) + kiểm "tổng bằng 0".

Độ phức tạp = **max(2¹²⁸·2¹³, 2⁶⁴·2¹³) = 2¹⁴¹** thao tác, dữ liệu **2¹³**. Con số **2¹⁴¹** là
**thực tế bất khả thi** → **ARADI đầy đủ KHÔNG bị phá**. Đây là kết quả *trung thực*, và việc nêu
rõ giới hạn là điểm mạnh, không phải điểm yếu.

## 5. Biến thể 6 vòng đã sửa key schedule — tấn công chạy được trong vài giây

Để **minh hoạ phương pháp chạy thực tế**, ta xét một **biến thể đã sửa lịch khoá (key schedule)**:

- Master key 256-bit được chiếu xuống `K_base = K₀ ⊕ … ⊕ K₇` (**32 bit**), tách thành
  `u` (16 bit cao) và `ℓ` (16 bit thấp). Khoá vòng sinh từ `u`/`ℓ` qua các phép quay.
- **Entropy hiệu dụng chỉ còn 32 bit** ⇒ không gian tìm kiếm sụp đổ.

Tấn công trên **6 vòng**:

- **Giai đoạn 1:** dò `u` (2¹⁶ ứng viên), nhận `u` nếu tổng-cube của X và Z **byte-wise bằng
  nhau** (AABB). Tỉ lệ báo nhầm ~**2⁻³²**/ứng viên.
- **Giai đoạn 2:** dò `ℓ` (2¹⁶ ứng viên), nhận nếu **tổng bốn từ = 0**. Tỉ lệ báo nhầm ~2⁻¹²⁸.

**Kết quả thực nghiệm** ([`../python/run_multicube_100.py`](../python/run_multicube_100.py)):

| Chỉ số | Giá trị |
|---|---|
| Master key ngẫu nhiên | **100/100** phục hồi đúng `(u, ℓ)` (bản luận văn gốc: 20/20) |
| Báo nhầm (Mode A, AABB đầy đủ) | **0** |
| Thời gian (Giai đoạn 2, Python+NumPy, 1 luồng) | **~6 giây/khoá** |

> **Quan trọng:** ta phục hồi **cặp 32-bit `(u, ℓ)`**, **KHÔNG** phải master key 256-bit. Phép
> chiếu `K_base = ⊕Kᵢ` **không khả nghịch** (mất 224 bit thông tin) ⇒ suy ngược master key là
> **bất khả về mặt thông tin**, không phải "khó".

## 6. Tối ưu multi-cube

Chạy **K** cube **độc lập** rồi chỉ giữ ứng viên qua **mọi** cube: tỉ lệ báo nhầm giai đoạn 1
giảm từ 2⁻³² xuống **2⁻³²ᴷ** (Mode A) / 2⁻¹⁶ᴷ (Mode B), đổi lại cần K lần nhiều bản rõ chọn hơn.
Đây là phần tối ưu đóng góp trong luận văn.

## 7. Tự chạy lại

```bash
make test-c      # cài đặt C khớp test-vector NSA  (exit 0 = PASS)
make reproduce   # tấn công 100 khoá ngẫu nhiên trên biến thể 6 vòng
make milp        # kiểm chứng AABB bằng MILP
```

## 8. Kết luận & giới hạn

- AABB là **tính chất cấu trúc** của hàm vòng ARADI (đã kiểm chứng độc lập bằng MILP + C).
- ARADI đầy đủ vẫn **an toàn thực tế** (~2¹⁴¹) — **không bị phá**.
- Trên **biến thể 6 vòng đã sửa**, tấn công cube phục hồi cặp 32-bit `(u, ℓ)` trong vài giây,
  100/100 — nhưng **không** chạm tới master key.
- Hướng phát triển: tổng quát hoá AABB sang **họ cipher Toffoli/χ**, và so sánh **neural
  distinguisher vs AABB** trên cùng testbed (xem `../generalize/`, `../ml/`).
