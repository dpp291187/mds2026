# Đưa mã nguồn và dữ liệu lên GitHub

Gói dùng để đưa lên GitHub là **`MDS_Reproducibility_GitHub.zip`**. Đây là kho mã độc lập, gồm Python, dữ liệu, hình, tài liệu hướng dẫn và thông tin nguồn tham khảo; không cần LaTeX. Bài báo, mã LaTeX và thư mục nộp tạp chí được giữ trong gói project riêng.

## 1. Giải nén và kiểm tra

Giải nén ZIP trước. Trong thư mục vừa giải nén phải có ngay `README.md`, `requirements.txt`, `code/`, `data/`, `figures/`, `provenance/` và các tệp thông tin kèm theo. ZIP không bọc thêm một thư mục cấp ngoài.

Có thể mở thư mục này bằng PyCharm và chạy:

```bash
python code/revision6_checks.py
```

Để tự xuất hình từ dữ liệu đã có:

```bash
python -m pip install -r requirements.txt
python code/make_figures.py
python code/make_supplementary_figures.py
```

Lệnh vẽ hình chỉ đọc dữ liệu; không chạy lại phép vét cạn. Hình chính được xuất PDF, PNG 600 dpi, SVG và EPS trong `figures/`; các bản `Fig1` và `Fig2` dạng PDF/EPS cũng được tạo tại thư mục gốc. Hai hình phụ được xuất PDF. README tiếng Anh trong kho giải thích đầy đủ cách chạy lại các kiểm chứng toán học theo đúng thứ tự.

## 2. Tạo kho và tải nội dung đã giải nén

Tạo một repository GitHub mới trong tài khoản của tác giả, chẳng hạn đặt tên `mds-integer-lifts`. Để người đọc và phản biện truy cập được, đặt kho ở chế độ công khai khi sẵn sàng công bố.

Có thể tải các tệp và thư mục **đã giải nén** qua giao diện thêm tệp của GitHub. Cần giữ nguyên các thư mục `code/`, `data/`, `figures/` và `provenance/`; kiểm tra lại danh sách tệp trước khi xác nhận commit. Không chỉ tải duy nhất tệp ZIP, vì khi đó README và mã nguồn không xuất hiện như một repository thông thường.

Nếu dùng Git trên máy tính, tạo repository GitHub trống, mở terminal tại thư mục đã giải nén và chạy các lệnh dưới đây sau khi thay URL bằng địa chỉ thật:

```bash
git init
git add .
git commit -m "Add reproducibility code, data and figures"
git branch -M main
git remote add origin https://github.com/USERNAME/REPOSITORY.git
git push -u origin main
```

`USERNAME/REPOSITORY` là chỗ giữ chỗ, không phải địa chỉ kho đã được tạo. Đăng nhập bằng phương thức GitHub đang dùng trên máy; không đưa mật khẩu hay token vào các tệp của project. Tệp `.gitignore` có sẵn sẽ loại môi trường Python, bộ nhớ đệm và thư mục build khỏi commit.

## 3. Thay đường dẫn trong bài báo

Sau khi kho đã công khai và mở được từ liên kết thật, thay `https://github.com/USERNAME/REPOSITORY` ở phần **Data availability** của bài báo bằng URL repository vừa tạo. Nếu kho thực tế dùng tên hoặc tài khoản khác, lấy trực tiếp URL trên trang repository để tránh gõ sai.

Kiểm tra rằng trang kho hiển thị README tiếng Anh, tải được dữ liệu trong `data/`, có `data/revision6_checks.json` và mở được hai mã vẽ hình. Đây là phần địa chỉ trực tuyến tác giả cần cập nhật; không điền một URL dự kiến như thể kho đã tồn tại.

## 4. Các lưu ý về hồ sơ

Tên và email trong README/CITATION được ghi là **Phuc-Phan Duong — duongphucphan@actvn.edu.vn**. Không có DOI hay giấy phép sử dụng tự đặt. Hai giấy phép đi kèm nguồn Plonky3 chỉ áp dụng theo thông báo của nguồn đó; chúng không tự động trở thành giấy phép cho toàn bộ mã của bài báo. Nếu muốn cho phép sử dụng lại mã nguồn theo một giấy phép cụ thể, tác giả có thể chọn giấy phép phù hợp khi công bố kho.

Nên giữ lại ZIP phát hành ban đầu để đối chiếu checksum. Chạy lại toàn bộ thực nghiệm có thể cập nhật các trường thời gian/môi trường trong báo cáo; chạy lại hình có thể đổi metadata của PDF/SVG. Những thay đổi đó không tự động đồng nghĩa với thay đổi kết quả toán học.
