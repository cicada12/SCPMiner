# A Tool for Drug Discovery

**Reference:**
Reddy, A. S., Reddy, P. K., Mondal, A., & Priyakumar, U. D. (2021). Mining subgraph coverage patterns from graph transactions. *International Journal of Data Science and Analytics*, 13(2), 105–121. [https://doi.org/10.1007/s41060-021-00292-y](https://doi.org/10.1007/s41060-021-00292-y)

## Project Overview

This tool is designed to generate subgraph coverage patterns from chemical compounds to aid in drug discovery. By analyzing subgraphs within the chemical compound structures, the tool helps identify key patterns and features that are crucial for drug development. These patterns can be leveraged to predict the properties of compounds and assist researchers in discovering new potential drugs.

## Installation

Follow the steps below to set up the project:

### 1. Clone the repository

```bash
git clone <repository_url>
cd <repository_name>
```

### 2. Run the frontend server

Navigate to the `frontend` folder and run the following command to start the frontend server:

```bash
cd frontend
npm run dev
```

### 3. Run the backend server

Navigate to the `backend` folder and run the following command to start the backend server:

```bash
cd backend
uvicorn main::app --reload
```
