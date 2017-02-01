function [ G ] = vesselsCalibre(G, segm)

    % compute the euclidean distance transform of the complement of the
    % segmentation to capture a raw approximation to the vessel diameter
    raw_diameters = bwdist(imcomplement(segm), 'euclidean')+1;
    
    displacement = 3; % to compute vessel normals
    sample_factor = 2; % to sample some normals
    
    figure, imshow(segm);
    hold on
    
    % compute vessel Calibre for each link
    for i = 1 : length(G.link)

        % if we have more than only 2 points to consider
        if (G.link(i).point  > 2 * displacement)
        
            % en pts guardo todas las coordenadas de los puntos del link
            pts = zeros(length(G.link(i).point),2);
            [pts(:,1), pts(:,2)] = ind2sub([G.w G.l], G.link(i).point);
            % en next_pts tengo las mismas coordenadas, pero desplazadas de
            % acuerdo al valor de displacement
            next_pts = circshift(pts',[0 -displacement])'; 
            % en vector_norm tengo el di?metro estimado, aproximado a
            % partir del m?ximo de la transformada de distancia en el
            % segmento
            vector_norm = ceil(max(raw_diameters(G.link(i).point)));
            
            % elimino las puntas, que me quedan cruzadas
            pts = pts(displacement:end-displacement,:);
            next_pts = next_pts(displacement:end-displacement,:);

            % muestreo el esqueleto cada 2 puntos
            v_sampled = next_pts(mod(1:length(pts),sample_factor)==1,:) - pts(mod(1:length(pts),sample_factor)==1,:);
            if (size(v_sampled,1) < 3)
                v_sampled = next_pts - pts;
            else
                pts = pts(mod(1:length(pts),sample_factor)==1,:);
                next_pts = next_pts(mod(1:length(next_pts),sample_factor)==1,:);
            end
            
            % calculo los vectores perpendiculares al esqueleto
            norm_ = ones(size(v_sampled,1),1) * vector_norm;
            u_sampled = zeros(size(v_sampled));
            u_sampled(:,1) = norm_ ./ sqrt(1 + (v_sampled(:,1) ./ v_sampled(:,2)).^2);
            u_sampled(:,2) = -1 * (v_sampled(:,1) ./ v_sampled(:,2)) .* u_sampled(:,1);

            % en caso de que haya alg?n delta_y=0, se corrige:
            u_sampled(v_sampled(:,1)==0, 2) = 0;
            u_sampled(v_sampled(:,2)==0, 2) = vector_norm;
            
            % calculo los puntos ortogonales al actual
            orthogonal_a = round(pts + u_sampled);
            orthogonal_b = round(pts - u_sampled);
            
            quiver(pts(:,2), pts(:,1), u_sampled(:,2), u_sampled(:,1));
            hold on
            quiver(pts(:,2), pts(:,1), -u_sampled(:,2), -u_sampled(:,1));
            hold on
            
            % identifico los puntos dentro del vaso
            G.link(i).vessel_width = zeros(size(orthogonal_a,1),1);
            for j = 1 : size(orthogonal_a,1)
                
                % obtengo las etiquetas en el vector normal
                [line] = improfile(segm, [orthogonal_a(j,2) orthogonal_b(j,2)], [orthogonal_a(j,1) orthogonal_b(j,1)], double(vector_norm*2+1));
                % el centerline est? localizado en el punto vector_norm+1.
                % usando esa informaci?n, puedo identificar d?nde empieza y
                % donde termina
                line_idx = find(line==0) - (vector_norm+1);
                try
                    neg = line_idx(line_idx<0);
                    if (isempty(neg)), a = vector_norm+1; else a = find(max(neg)==line_idx); end
                    pos = line_idx(line_idx>0);
                    if (isempty(pos)), b = vector_norm+1; else b = find(min(pos)==line_idx); end
                catch ex
                    disp('A');
                end
                G.link(i).vessel_width(j) = abs(a-b)+1;
                %length(line)
                
            end
            mean(G.link(i).vessel_width)
            % calculo el ancho como la distancia euclidea entre el primer
            % punto marcado con 1 en cada recta
%             for j = 1 : size(orthogonal_b,1)
%                 
%                 % identify pixels inside the vessel
%                 [cx cy line] = improfile(segm, [orthogonal_a(j,1) orthogonal_b(j,1)], [orthogonal_a(j,2) orthogonal_b(j,2)]);
%                 
%                 % identify the closest edge to the centerline
%                 vessel_pixels = (find(line==1) - 51);
%                 first_piece = find(vessel_pixels < 0);
%                 if (~isempty(first_piece))
%                     [~, idx_begin] = max(first_piece);
%                 else
%                     idx_begin = 51;
%                 end
%                 last_piece = find(vessel_pixels < 0);
%                 if (~isempty(last_piece))
%                     [~, idx_end] = min(last_piece);
%                 else
%                     idx_end = 51;
%                 end
%                 
%                 
%                 
%                 
%                 line = find(1-segm(sub2ind([G.w G.l], x, y)));
%                 first_point = ind2sub([G.w G.l], line(1));
%                 
%                 [x, y]=bresenham(pts(j,1),pts(j,2),orthogonal_a(j,1),orthogonal_a(j,2));
%                 line = find(1-segm(sub2ind([G.w G.l], x, y)));
%                 end_point = ind2sub([G.w G.l], line(1));
%                 
%             end     
            
        end   
    end

% 
% 

end    
    
    



