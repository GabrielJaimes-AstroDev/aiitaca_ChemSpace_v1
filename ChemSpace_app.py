import streamlit as st
import pandas as pd
import numpy as np
import plotly.express as px
import plotly.graph_objects as go
from plotly.subplots import make_subplots
import json
import re

# Page configuration
st.set_page_config(
    page_title="Chemical Space Visualization",
    page_icon="🧪",
    layout="wide",
    initial_sidebar_state="expanded"
)

# Custom CSS
st.markdown("""
<style>
    .main-header {
        font-size: 2.5rem;
        color: #1f77b4;
        text-align: center;
        margin-bottom: 2rem;
    }
    .sub-header {
        font-size: 1.5rem;
        color: #2ca02c;
        margin-bottom: 1rem;
    }
    .info-box {
        background-color: #f0f8ff;
        padding: 1rem;
        border-radius: 0.5rem;
        border-left: 4px solid #1f77b4;
        margin-bottom: 1rem;
    }
    .molecule-card {
        background-color: #f8f9fa;
        padding: 0.5rem;
        border-radius: 0.3rem;
        margin: 0.2rem 0;
        border-left: 3px solid #9467bd;
    }
    .priority-db { border-left-color: #ff7f0e; }
    .candidate { border-left-color: #ff6b6b; }
    .selected { border-left-color: #dc3545; }
</style>
""", unsafe_allow_html=True)

def load_plot_data(file_path):
    """Load plot data from JSON file"""
    try:
        with open(file_path, 'r') as f:
            return json.load(f)
    except Exception as e:
        st.error(f"Error loading plot data: {str(e)}")
        return None

def format_chemical_formula(formula):
    """Format chemical formula with subscripts"""
    if not formula or pd.isna(formula):
        return ""
    
    formula = re.sub(r'(\d+)', r'<sub>\1</sub>', str(formula))
    formula = re.sub(r'([+-]\d*)', r'<sup>\1</sup>', formula)
    
    return formula

def create_interactive_plot(plot_data):
    """Create interactive Plotly visualization with proper layer ordering"""
    
    # Create DataFrame from plot data
    df_points = pd.DataFrame({
        'x': [point[0] for point in plot_data['points']],
        'y': [point[1] for point in plot_data['points']],
        'label': plot_data['labels'],
        'database': plot_data['databases'],
        'formula': plot_data['formulas'],
        'color': plot_data['colors']
    })
    
    # Create the main scatter plot
    fig = go.Figure()
    
    # Define z-order for proper layering
    zorders = {
        'PubChem': 1,
        'KIDA': 2,
        'All_Discoveries': 3,
        'TMC_1': 4,
        'GUAPOS': 5,
        'candidate': 6,
        'selected': 7,
        'guapos': 8,
        'connections': 3
    }
    
    # Add PubChem first (background)
    pubchem_mask = df_points['database'] == 'PubChem'
    if pubchem_mask.any():
        pubchem_data = df_points[pubchem_mask]
        hover_text = []
        for _, row in pubchem_data.iterrows():
            text = f"Database: {row['database']}<br>"
            if row['label']:
                text += f"Name: {row['label']}<br>"
            if row['formula']:
                text += f"Formula: {format_chemical_formula(row['formula'])}"
            hover_text.append(text)
        
        fig.add_trace(go.Scatter(
            x=pubchem_data['x'],
            y=pubchem_data['y'],
            mode='markers',
            marker=dict(
                color=pubchem_data['color'].iloc[0],
                size=6,
                opacity=0.2
            ),
            name='PubChem',
            hovertext=hover_text,
            hoverinfo='text',
            showlegend=True,
            legendgroup='databases'
        ))
    
    # Add other databases
    other_dbs = [db for db in df_points['database'].unique() if db != 'PubChem']
    for db in other_dbs:
        db_mask = df_points['database'] == db
        db_data = df_points[db_mask]
        
        hover_text = []
        for _, row in db_data.iterrows():
            text = f"Database: {row['database']}<br>"
            if row['label']:
                text += f"Name: {row['label']}<br>"
            if row['formula']:
                text += f"Formula: {format_chemical_formula(row['formula'])}"
            hover_text.append(text)
        
        fig.add_trace(go.Scatter(
            x=db_data['x'],
            y=db_data['y'],
            mode='markers',
            marker=dict(
                color=db_data['color'].iloc[0],
                size=10,
                opacity=0.7
            ),
            name=db,
            hovertext=hover_text,
            hoverinfo='text',
            showlegend=True,
            legendgroup='databases'
        ))
    
    # Add candidate neighbors if they exist
    if plot_data.get('candidate_neighbors'):
        candidate_indices = set()
        for candidate in plot_data['candidate_neighbors']:
            if not candidate['is_priority']:
                candidate_indices.add(candidate['neighbor_index'])
        
        candidate_indices = list(candidate_indices)
        if candidate_indices:
            candidate_points = []
            candidate_info = []
            for idx in candidate_indices:
                if idx < len(plot_data['points']):
                    candidate_points.append(plot_data['points'][idx])
                    candidate_info.append({
                        'database': 'Candidate',
                        'formula': plot_data['formulas'][idx],
                        'label': plot_data['labels'][idx]
                    })
            
            if candidate_points:
                hover_text = []
                for info in candidate_info:
                    text = "Type: Candidate Neighbor<br>"
                    if info['label']:
                        text += f"Name: {info['label']}<br>"
                    if info['formula']:
                        text += f"Formula: {format_chemical_formula(info['formula'])}"
                    hover_text.append(text)
                
                fig.add_trace(go.Scatter(
                    x=[p[0] for p in candidate_points],
                    y=[p[1] for p in candidate_points],
                    mode='markers',
                    marker=dict(
                        color='lightcoral',
                        size=8,
                        opacity=0.6
                    ),
                    name='Candidate Neighbors',
                    hovertext=hover_text,
                    hoverinfo='text',
                    showlegend=True,
                    legendgroup='neighbors'
                ))
    
    # Add selected neighbors if they exist
    if plot_data.get('knn_connections'):
        selected_indices = set()
        for connection in plot_data['knn_connections']:
            if connection['type'] == 'selected':
                selected_indices.add(connection['neighbor_index'])
        
        selected_indices = list(selected_indices)
        if selected_indices:
            selected_points = []
            selected_info = []
            for idx in selected_indices:
                if idx < len(plot_data['points']):
                    selected_points.append(plot_data['points'][idx])
                    selected_info.append({
                        'database': plot_data['databases'][idx],
                        'formula': plot_data['formulas'][idx],
                        'label': plot_data['labels'][idx]
                    })
            
            if selected_points:
                hover_text = []
                for info in selected_info:
                    text = "Type: Selected Neighbor<br>"
                    text += f"Database: {info['database']}<br>"
                    if info['label']:
                        text += f"Name: {info['label']}<br>"
                    if info['formula']:
                        text += f"Formula: {format_chemical_formula(info['formula'])}"
                    hover_text.append(text)
                
                fig.add_trace(go.Scatter(
                    x=[p[0] for p in selected_points],
                    y=[p[1] for p in selected_points],
                    mode='markers',
                    marker=dict(
                        color='red',
                        size=10,
                        opacity=0.8
                    ),
                    name='Selected Neighbors',
                    hovertext=hover_text,
                    hoverinfo='text',
                    showlegend=True,
                    legendgroup='neighbors'
                ))
    
    # Add GUAPOS molecules if they exist
    guapos_indices = [i for i, db in enumerate(plot_data['databases']) if db == 'GUAPOS']
    if guapos_indices:
        guapos_points = [plot_data['points'][i] for i in guapos_indices]
        guapos_info = [{
            'formula': plot_data['formulas'][i],
            'label': plot_data['labels'][i]
        } for i in guapos_indices]
        
        hover_text = []
        for info in guapos_info:
            text = "Database: GUAPOS<br>"
            if info['label']:
                text += f"Name: {info['label']}<br>"
            if info['formula']:
                text += f"Formula: {format_chemical_formula(info['formula'])}"
            hover_text.append(text)
        
        fig.add_trace(go.Scatter(
            x=[p[0] for p in guapos_points],
            y=[p[1] for p in guapos_points],
            mode='markers',
            marker=dict(
                color='blue',
                size=12,
                opacity=0.9,
                line=dict(color='black', width=1)
            ),
            name='GUAPOS Molecules',
            hovertext=hover_text,
            hoverinfo='text',
            showlegend=True,
            legendgroup='guapos'
        ))
    
    # Add KNN connections if they exist
    if plot_data.get('knn_connections'):
        for connection in plot_data['knn_connections']:
            if (connection['guapos_index'] < len(plot_data['points']) and 
                connection['neighbor_index'] < len(plot_data['points'])):
                
                x0, y0 = plot_data['points'][connection['guapos_index']]
                x1, y1 = plot_data['points'][connection['neighbor_index']]
                
                line_style = dict(
                    color='red' if connection['type'] == 'selected' else 'gray',
                    width=2 if connection['type'] == 'selected' else 1,
                    dash='solid' if connection['type'] == 'selected' else 'dash'
                )
                
                fig.add_trace(go.Scatter(
                    x=[x0, x1],
                    y=[y0, y1],
                    mode='lines',
                    line=line_style,
                    opacity=0.6,
                    showlegend=False,
                    hoverinfo='skip'
                ))
    
    # Update layout
    fig.update_layout(
        title="Chemical Space Visualization - PCA/UMAP Projection",
        xaxis_title="Dimension 1",
        yaxis_title="Dimension 2",
        hovermode='closest',
        showlegend=True,
        legend=dict(
            orientation="h",
            yanchor="bottom",
            y=1.02,
            xanchor="right",
            x=1
        ),
        height=700,
        width=1000
    )
    
    fig.update_xaxes(showgrid=False, zeroline=False, showticklabels=False)
    fig.update_yaxes(showgrid=False, zeroline=False, showticklabels=False)
    
    return fig

def create_database_summary(plot_data):
    """Create database summary statistics"""
    databases = plot_data['databases']
    db_counts = pd.Series(databases).value_counts().reset_index()
    db_counts.columns = ['Database', 'Count']
    
    fig = px.bar(
        db_counts,
        x='Database',
        y='Count',
        color='Database',
        title="Molecule Distribution by Database",
        color_discrete_map={
            'GUAPOS': '#1f77b4',
            'TMC_1': '#ff7f0e',
            'All_Discoveries': '#2ca02c',
            'KIDA': '#FFFF00',
            'PubChem': '#9467bd'
        }
    )
    
    fig.update_layout(
        xaxis_title="Database",
        yaxis_title="Number of Molecules",
        showlegend=False
    )
    
    return fig

def create_knn_summary(plot_data):
    """Create KNN analysis summary"""
    if not plot_data.get('knn_connections'):
        return None
    
    connection_types = [conn['type'] for conn in plot_data['knn_connections']]
    type_counts = pd.Series(connection_types).value_counts().reset_index()
    type_counts.columns = ['Type', 'Count']
    
    fig = px.pie(
        type_counts,
        values='Count',
        names='Type',
        title="KNN Connection Types",
        color='Type',
        color_discrete_map={
            'selected': 'red',
            'candidate': 'lightcoral'
        }
    )
    
    return fig

def main():
    st.markdown('<h1 class="main-header">🧪 Chemical Space Visualization</h1>', unsafe_allow_html=True)
    
    # Sidebar
    with st.sidebar:
        st.header("📊 Data Configuration")
        
        uploaded_file = st.file_uploader(
            "Upload plot data JSON",
            type=['json'],
            help="Upload the plot_data.json file generated by the analysis script"
        )
        
        if uploaded_file is not None:
            try:
                plot_data = json.load(uploaded_file)
                st.success("✅ Data loaded successfully!")
            except Exception as e:
                st.error(f"❌ Error loading file: {str(e)}")
                plot_data = None
        else:
            try:
                plot_data = load_plot_data("plot_data.json")
                if plot_data:
                    st.info("📁 Using default plot_data.json file")
            except:
                plot_data = None
                st.warning("⚠️ Please upload a plot data JSON file")
        
        if plot_data:
            st.markdown("---")
            st.header("🎨 Display Options")
            
            # Database summary
            db_counts = pd.Series(plot_data['databases']).value_counts()
            st.write("**Database Summary:**")
            for db, count in db_counts.items():
                st.write(f"• {db}: {count} molecules")
            
            # KNN info
            if plot_data.get('knn_connections'):
                selected = sum(1 for conn in plot_data['knn_connections'] if conn['type'] == 'selected')
                candidate = sum(1 for conn in plot_data['knn_connections'] if conn['type'] == 'candidate')
                st.write(f"• KNN Connections: {len(plot_data['knn_connections'])}")
                st.write(f"  - Selected: {selected}")
                st.write(f"  - Candidate: {candidate}")
            
            st.markdown("---")
            st.header("🔍 Interaction Guide")
            st.info("""
            - **Hover** over points to see molecule details
            - **Click** on points to view detailed information
            - **Zoom** with mouse wheel or touchpad
            - **Pan** by dragging the plot
            - **Reset view** with home button in toolbar
            - **PubChem** molecules are in the background
            - **Priority databases** are shown on top
            """)
    
    # Main content
    if plot_data is None:
        st.warning("""
        ## Welcome to Chemical Space Visualization!
        
        To get started:
        
        1. Run the main analysis script to generate `plot_data.json`
        2. Upload the JSON file using the sidebar
        3. Or place `plot_data.json` in the same directory as this app
        
        The enhanced visualization now includes:
        - Candidate neighbors (light red)
        - Selected neighbors (red)
        - Proper layer ordering with PubChem in background
        - KNN connection types
        """)
        return
    
    # Create tabs
    tab1, tab2, tab3, tab4 = st.tabs(["📈 Interactive Map", "📊 Statistics", "🔍 Molecule Explorer", "📋 KNN Analysis"])
    
    with tab1:
        st.markdown('<h2 class="sub-header">Interactive Chemical Space Map</h2>', unsafe_allow_html=True)
        
        fig = create_interactive_plot(plot_data)
        st.plotly_chart(fig, use_container_width=True)
        
        col1, col2 = st.columns(2)
        
        with col1:
            st.markdown("""
            <div class="info-box">
            <h4>Enhanced Visualization</h4>
            <p>• <strong style="color:lightcoral;">Candidate neighbors</strong>: All potential matches<br>
            • <strong style="color:red;">Selected neighbors</strong>: Final choices after filtering<br>
            • <strong>PubChem background</strong>: Proper layer ordering<br>
            • <strong>Solid red lines</strong>: Selected connections<br>
            • <strong>Dashed gray lines</strong>: Candidate connections</p>
            </div>
            """, unsafe_allow_html=True)
        
        with col2:
            st.markdown("""
            <div class="info-box">
            <h4>Layer Ordering</h4>
            <p>Molecules are layered from bottom to top:<br>
            1. PubChem (background)<br>
            2. Other databases<br>
            3. Candidate neighbors<br>
            4. Selected neighbors<br>
            5. GUAPOS molecules<br>
            6. Connection lines</p>
            </div>
            """, unsafe_allow_html=True)
    
    with tab2:
        st.markdown('<h2 class="sub-header">Database Statistics</h2>', unsafe_allow_html=True)
        
        col1, col2 = st.columns(2)
        
        with col1:
            db_fig = create_database_summary(plot_data)
            st.plotly_chart(db_fig, use_container_width=True)
        
        with col2:
            total_molecules = len(plot_data['points'])
            unique_databases = len(set(plot_data['databases']))
            
            st.metric("Total Molecules", f"{total_molecules:,}")
            st.metric("Unique Databases", unique_databases)
            
            if plot_data.get('knn_connections'):
                selected = sum(1 for conn in plot_data['knn_connections'] if conn['type'] == 'selected')
                candidate = sum(1 for conn in plot_data['knn_connections'] if conn['type'] == 'candidate')
                st.metric("Total KNN Connections", len(plot_data['knn_connections']))
                st.metric("Selected Connections", selected)
                st.metric("Candidate Connections", candidate)
    
    with tab3:
        st.markdown('<h2 class="sub-header">Molecule Explorer</h2>', unsafe_allow_html=True)
        
        molecules_df = pd.DataFrame({
            'Database': plot_data['databases'],
            'Name': plot_data['labels'],
            'Formula': plot_data['formulas'],
            'X': [p[0] for p in plot_data['points']],
            'Y': [p[1] for p in plot_data['points']]
        })
        
        # Search and filter
        col1, col2 = st.columns(2)
        
        with col1:
            selected_db = st.multiselect(
                "Filter by Database",
                options=sorted(molecules_df['Database'].unique()),
                default=sorted(molecules_df['Database'].unique())
            )
        
        with col2:
            search_term = st.text_input("Search by Name or Formula", "")
        
        # Apply filters
        filtered_df = molecules_df[molecules_df['Database'].isin(selected_db)]
        if search_term:
            filtered_df = filtered_df[
                filtered_df['Name'].str.contains(search_term, case=False, na=False) |
                filtered_df['Formula'].str.contains(search_term, case=False, na=False)
            ]
        
        st.write(f"**Found {len(filtered_df)} molecules**")
        
        # Pagination
        page_size = 20
        total_pages = max(1, (len(filtered_df) + page_size - 1) // page_size)
        page = st.number_input("Page", min_value=1, max_value=total_pages, value=1)
        
        start_idx = (page - 1) * page_size
        end_idx = min(start_idx + page_size, len(filtered_df))
        
        # Display molecules
        for idx in range(start_idx, end_idx):
            molecule = filtered_df.iloc[idx]
            with st.expander(f"{molecule['Name'] or 'Unknown'} - {molecule['Database']}"):
                col1, col2 = st.columns(2)
                
                with col1:
                    st.write(f"**Database:** {molecule['Database']}")
                    st.write(f"**Name:** {molecule['Name'] or 'Unknown'}")
                    st.write(f"**Formula:** {format_chemical_formula(molecule['Formula'])}", unsafe_allow_html=True)
                
                with col2:
                    st.write(f"**Coordinates:** ({molecule['X']:.3f}, {molecule['Y']:.3f})")
                    st.write(f"**Index:** {filtered_df.index[idx]}")
        
        # Download option
        csv = filtered_df.to_csv(index=False)
        st.download_button(
            label="📥 Download Filtered Data",
            data=csv,
            file_name="filtered_molecules.csv",
            mime="text/csv"
        )
    
    with tab4:
        st.markdown('<h2 class="sub-header">KNN Analysis Details</h2>', unsafe_allow_html=True)
        
        if not plot_data.get('candidate_neighbors'):
            st.info("No KNN analysis data available in this file.")
            return
        
        # KNN summary
        col1, col2 = st.columns(2)
        
        with col1:
            knn_fig = create_knn_summary(plot_data)
            if knn_fig:
                st.plotly_chart(knn_fig, use_container_width=True)
        
        with col2:
            candidate_data = plot_data['candidate_neighbors']
            total_candidates = len(candidate_data)
            priority_candidates = sum(1 for c in candidate_data if c['is_priority'])
            valid_candidates = total_candidates - priority_candidates
            
            st.metric("Total Candidate Pairs", total_candidates)
            st.metric("Priority Database Candidates", priority_candidates)
            st.metric("Valid Candidates", valid_candidates)
            
            # Show priority databases
            if plot_data.get('priority_databases'):
                st.write("**Priority Databases (excluded):**")
                for db in plot_data['priority_databases']:
                    st.write(f"• {db}")
        
        # Show candidate details
        st.write("**Candidate Neighbor Details:**")
        
        # Group by GUAPOS molecule
        guapos_candidates = {}
        for candidate in candidate_data:
            guapos_idx = candidate['guapos_index']
            if guapos_idx not in guapos_candidates:
                guapos_candidates[guapos_idx] = []
            guapos_candidates[guapos_idx].append(candidate)
        
        for guapos_idx, candidates in list(guapos_candidates.items())[:5]:  # Show first 5
            if guapos_idx < len(plot_data['points']):
                guapos_name = plot_data['labels'][guapos_idx] or f"GUAPOS_{guapos_idx}"
                with st.expander(f"GUAPOS Molecule: {guapos_name}"):
                    for candidate in sorted(candidates, key=lambda x: x['rank']):
                        status = "PRIORITY" if candidate['is_priority'] else "CANDIDATE"
                        status_color = "priority-db" if candidate['is_priority'] else "candidate"
                        if candidate['type'] == 'selected':
                            status_color = "selected"
                        
                        neighbor_idx = candidate['neighbor_index']
                        if neighbor_idx < len(plot_data['points']):
                            neighbor_name = plot_data['labels'][neighbor_idx] or f"Neighbor_{neighbor_idx}"
                            neighbor_db = plot_data['databases'][neighbor_idx]
                            
                            st.markdown(f"""
                            <div class="molecule-card {status_color}">
                                <strong>Rank {candidate['rank']}:</strong> {neighbor_name}<br>
                                <strong>Database:</strong> {neighbor_db}<br>
                                <strong>Status:</strong> {status}<br>
                                <strong>Distance:</strong> {candidate['distance']:.4f}
                            </div>
                            """, unsafe_allow_html=True)

if __name__ == "__main__":
    main()
