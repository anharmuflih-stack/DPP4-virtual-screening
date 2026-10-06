import os

filepath = '/Users/aanmuf/drug_discovery/dpp4_project/scripts/11_dashboard_ui.py'
with open(filepath, 'r') as f:
    content = f.read()

css = """
st.markdown('''
<style>
/* Modern Gradient Background */
.stApp {
    background: linear-gradient(135deg, #0f172a 0%, #1e1b4b 100%);
    color: #e2e8f0;
}

/* Glassmorphism Cards */
div.css-1r6slb0, div.css-12oz5g7 {
    background: rgba(255, 255, 255, 0.05);
    backdrop-filter: blur(10px);
    border-radius: 15px;
    border: 1px solid rgba(255, 255, 255, 0.1);
    padding: 20px;
    box-shadow: 0 8px 32px 0 rgba(0, 0, 0, 0.3);
}

/* Premium Buttons */
.stButton>button {
    background: linear-gradient(90deg, #3b82f6 0%, #8b5cf6 100%);
    color: white;
    border: none;
    border-radius: 8px;
    padding: 10px 24px;
    font-weight: 600;
    transition: all 0.3s ease;
    box-shadow: 0 4px 15px rgba(139, 92, 246, 0.3);
}
.stButton>button:hover {
    transform: translateY(-2px);
    box-shadow: 0 6px 20px rgba(139, 92, 246, 0.5);
    border: none;
    color: white;
}

/* Typography Enhancements */
h1, h2, h3 {
    color: #f8fafc !important;
    font-family: 'Inter', sans-serif;
    font-weight: 700;
    letter-spacing: -0.5px;
}

/* Text Inputs */
.stTextInput>div>div>input {
    background-color: rgba(0,0,0,0.2) !important;
    color: #60a5fa !important;
    border: 1px solid rgba(255,255,255,0.1) !important;
    border-radius: 8px;
    font-family: 'Courier New', monospace;
    font-weight: bold;
}

/* Metric text */
.stMarkdown p {
    font-size: 1.05rem;
}

/* Hide default streamlit menu */
#MainMenu {visibility: hidden;}
footer {visibility: hidden;}
</style>
''', unsafe_allow_html=True)
"""

if "st.markdown('''" not in content:
    # Insert right after st.set_page_config
    parts = content.split('st.title')
    new_content = parts[0] + css + '\nst.title' + parts[1]
    with open(filepath, 'w') as f:
        f.write(new_content)
        
print("UX styles added.")
