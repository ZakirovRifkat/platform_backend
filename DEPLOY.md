# Deploy

## 1) Подготовка сервера (Ubuntu)

```bash
sudo apt update
sudo apt install -y python3 python3-venv python3-pip
```

## 2) Клонирование и установка зависимостей

```bash
cd /opt
sudo git clone <YOUR_REPO_URL> platform_backend
cd platform_backend
python3 -m venv .venv
source .venv/bin/activate
pip install --upgrade pip
pip install -r requirements.txt
```

## 3) Проверка ручного запуска

```bash
source /opt/platform_backend/.venv/bin/activate
uvicorn main:app --host 0.0.0.0 --port 8000
```

Проверка:
- `http://SERVER_IP:8000/health`
- `http://SERVER_IP:8000/docs`

## 4) Автозапуск через systemd (без Nginx)

Создай файл `/etc/systemd/system/platform-backend.service`:

```ini
[Unit]
Description=Platform Backend FastAPI
After=network.target

[Service]
User=root
Group=root
WorkingDirectory=/opt/platform_backend
Environment="PYTHONUNBUFFERED=1"
ExecStart=/opt/platform_backend/.venv/bin/uvicorn main:app --host 0.0.0.0 --port 8000 --workers 2
Restart=always
RestartSec=3

[Install]
WantedBy=multi-user.target
```

Применить:

```bash
sudo systemctl daemon-reload
sudo systemctl enable platform-backend
sudo systemctl start platform-backend
sudo systemctl status platform-backend
```

Логи:

```bash
journalctl -u platform-backend -f
```

## 5) Открыть порт в firewall

Если используется UFW:

```bash
sudo ufw allow 8000/tcp
sudo ufw reload
sudo ufw status
```

## 6) Обновления

```bash
cd /opt/platform_backend
git pull
source .venv/bin/activate
pip install -r requirements.txt
sudo systemctl restart platform-backend
```

## Пример запроса

```bash
curl "http://SERVER_IP:8000/model?A=0.5&B=0.2&C=0.3&delta=0.1&initial=1&initial=0&initial=0&initial=0"
```
