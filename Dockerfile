FROM python:3.11-slim

# Install system dependencies
RUN apt-get update && apt-get install -y curl && rm -rf /var/lib/apt/lists/*

# Set working directory
WORKDIR /app

# Install Python dependencies
COPY requirements.txt .
RUN pip install --no-cache-dir -r requirements.txt

# Copy the rest of the project
COPY . .

# Download Vina for Linux
RUN curl -L https://github.com/ccsb-scripps/AutoDock-Vina/releases/download/v1.2.5/vina_1.2.5_linux_x86_64 -o scripts/vina_linux
RUN chmod +x scripts/vina_linux

# Make directories writable for dynamic outputs during docking
RUN mkdir -p /app/results
RUN chmod -R 777 /app/results /app/data

# Expose port 7860 (Hugging Face Spaces default)
EXPOSE 7860

# Run with Gunicorn on port 7860
CMD ["gunicorn", "scripts.12_flask_app:app", "--bind", "0.0.0.0:7860", "--timeout", "120"]
