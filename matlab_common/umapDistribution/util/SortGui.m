%   AUTHORSHIP
%   Primary Developer: Stephen Meehan <swmeehan@stanford.edu> 
%   Math Lead & Secondary Developer:  Connor Meehan <cgmeehan@alumni.caltech.edu>
%   Bioinformatics Lead:  Wayne Moore <wmoore@stanford.edu>
%   Provided by the Herzenberg Lab at Stanford University 
%   License: BSD 3 clause
%
classdef SortGui < handle
    properties(Constant)
        PROP_ORDER='sortGui.sortOrder';
        PROP_ASCENDING='sortGui.sortAscending';
    end
    
    properties
        fncRefresh;
        dlg;
    end
    
    properties(SetAccess=private)
        path;
        strings;
        allChb;
        allMsg;
        sortIdxs;
        sortTypes;
        sortKeys;
        btnSortDir;
        btnSearch;
        options;
        sortAscending=true;
        nOptions;
        checkBoxes
        innerPnl;
        cb;
        map;
        propPrefix;
        allChbPnl;
        searchPnl;
        searchJtf;
        searching=false;
    end
    
    methods
        function this=SortGui(dlg, allChb, allMsg, allChbPnl, options, ...
                checkBoxes, innerPnl)
            pnl=Gui.Panel;
            this.dlg=dlg;
            this.allChb=allChb;
            this.allMsg=allMsg;
            this.options=options;
            this.nOptions=length(options);
            this.sortIdxs=1:this.nOptions;
            this.checkBoxes=checkBoxes;
            this.innerPnl=innerPnl;
            [this.sortKeys, this.sortTypes]=Html.DecodeSort(options{1});
            if ~isempty(this.sortKeys)
                jl=Gui.Label('Sort:');
                pnl.add(jl);
                fl=javaObjectEDT('java.awt.FlowLayout', 0,0,0);
                pnl.setLayout(fl);
                labels=['original' this.sortKeys];
                cb=Gui.Combo(labels, 1,[],[],@(h,e)sortItems(this, h), ...
                    'Change order');
                pnl.add(cb);
                this.cb=cb;
                allChbPnl.add(Gui.Label('  '), 'Center');
            end
            this.path=BasicMap.Path;
            this.btnSortDir=Gui.ImageButton(...
                fullfile(this.path, 'sortDown16.png'), ...
                'Order is ascending, click for DESCENDING?', ...
                @(h,e)toggleAscending(this));
            pnl.add(this.btnSortDir);
            this.btnSearch=Gui.ImageButton(...
                fullfile(this.path, 'magnify.png'), ...
                'Change order to descending', ...
                @(h,e)toggleSearch(this));
            pnl.add(this.btnSearch);
            this.initSearch;
            allChbPnl.add(pnl, 'East');
            this.allChbPnl=allChbPnl;
        end
        
        function setSortDir(this)
            if this.sortAscending
                this.btnSortDir.setIcon(Gui.Icon('sortDown16.png'));
                this.btnSortDir.setToolTipText(...
                    'Order is ascending, click for DESCENDING.');
            else
                this.btnSortDir.setIcon(Gui.Icon('sortUp16.png'));
                this.btnSortDir.setToolTipText(...
                    'Order is descending, click for ASCENDING.');
            end
        end
        
        function toggleAscending(this)
            this.sortIdxs=flip(this.sortIdxs);
            if ~isempty(this.fncRefresh)
                feval(this.fncRefresh);
            else
                this.innerPnl.removeAll;
                for jj=1:this.nOptions
                    idx=this.sortIdxs(jj);
                    this.innerPnl.add(this.checkBoxes.get(idx-1));
                end
            end
            this.sortAscending=~this.sortAscending;
            if ~isempty(this.map)
                this.map.set(this.propAscending, this.sortAscending);
            end
            this.setSortDir;
            if ~isempty(this.dlg)
                this.dlg.pack;
            end
        end
        
        function prop=propOrder(this)
            prop=[this.propPrefix SortGui.PROP_ORDER];
        end
        
        function prop=propAscending(this)
            prop=[this.propPrefix SortGui.PROP_ASCENDING];
        end
        
        function setProperties(this, map, propPrefix)
            this.map=map;
            this.propPrefix=propPrefix;
            if map.has(this.propOrder)
                idx=map.getNumeric(this.propOrder);
                this.cb.setSelectedIndex(idx-1);
                %this.sort(idx);
            end
            if ~map.is(this.propAscending, true)
                this.toggleAscending;
                this.setSortDir;
            end
        end
        
        function sortItems(this, h)
            idx=h.getSelectedIndex+1;
            if ~isempty(this.map)
                this.map.set(this.propOrder, idx);
            end
            this.sort(idx);
        end
        
        function sort(this, idx)
            if idx==1
                this.sortIdxs=1:this.nOptions;
            else
                isNum=isequal('N', this.sortTypes{idx-1});
                if isNum
                    values=zeros(1, this.nOptions);
                else
                    values=cell(1, this.nOptions);
                end
                key=this.sortKeys{idx-1};
                for jj=1:this.nOptions
                    v=Html.DecodeSortValues(this.options{jj}, key);
                    if isNum
                        values(jj)=str2double(v{1});
                    else
                        if isempty(v)
                            values{jj}='';
                        else
                            values{jj}=v{1};
                        end
                    end
                end
                [~, this.sortIdxs]=sort(values);
            end
            if this.sortAscending
                this.sortIdxs=flip(this.sortIdxs);
            end
            this.sortAscending=~this.sortAscending;
            this.toggleAscending;
            if ~isempty(this.map)
                this.map.set(this.propAscending, this.sortAscending);
            end
        end
        
        function initSearch(this)
            jtf=javaObjectEDT('javax.swing.JTextField');
            jtf.setColumns(15)
            %jtf.setHorizontalAlignment(jtf.RIGHT);
            jj=handle(jtf, 'CallbackProperties');
            set(jj, 'KeyTypedCallback', @(h,e)search(this, h));
            this.searchPnl=Gui.BorderPanel;
            pnl=Gui.Panel;
            pnl.add(jtf);
            btnClear=Gui.ImageButton(...
                fullfile(BasicMap.Path, 'cancel.gif'), ...
                'Clear search', @(h,e)clearSearch(this));
            pnl.add(btnClear);
            this.searchPnl.add(pnl, 'East');
            this.searchJtf=jtf;
        end
        
        function toggleSearch(this)
            if this.searching
                this.allChbPnl.remove(this.searchPnl);
            else
                this.allChbPnl.add(this.searchPnl, 'South');
                this.searchJtf.requestFocus;
            end
            this.searching=~this.searching;
            this.dlg.pack;
        end
        
        function clearSearch(this)
            this.searchJtf.setText('');
            this.search(this.searchJtf);
        end
        
        function search(this, h)
            N=this.nOptions-1;
            s=lower(char(h.getText));
            N=this.nOptions-1;
            if isempty(this.strings)
                this.strings={};
                for i=0:N
                    this.strings{end+1}=string(lower(...
                        char(edu.stanford.facs.swing.Basics.RemoveXml(...
                        this.checkBoxes.get(i).getText))));
                end
            end
            if isempty(s)
                for i=0:N
                    cb_=this.checkBoxes.get(i);
                    cb_.setEnabled(true);
                end
            else
                first=false;
                N=this.nOptions-1;
                for i=0:N
                    cb_=this.checkBoxes.get(i);
                    if this.strings{i+1}.contains(s)
                        cb_.setEnabled(true);
                        if ~first
                            cb_.scrollRectToVisible(cb_.getBounds)
                            first=true;
                        end
                    else
                        cb_.setEnabled(false);
                    end
                end
            end
            %this.dlg.pack;
        end
        
        function [idxs_, N]=setAllChbText(this)
            [idxs_, N]=Gui.SetAllChb(this.allChb, this.allMsg, this.checkBoxes);
        end
    end
end


