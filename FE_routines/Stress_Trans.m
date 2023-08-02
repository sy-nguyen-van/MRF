function [pred_x_Unit18,pred_y_Unit18,pred_xy_Unit18] = Stress_Trans(pred_x, pred_y, pred_xy)
global  TPMS SGs FE
switch TPMS
    case 'Primitive'
        %---------U1-U3-U6-U8
        pred_x_U1 = pred_x(1:SGs.ele_SGs/2,:,:);
        pred_y_U1 = pred_y(1:SGs.ele_SGs/2,:,:);
        pred_xy_U1 = pred_xy(1:SGs.ele_SGs/2,:,:);
        
        pred_x_U3 = pred_x_U1;
        pred_y_U3 = pred_y_U1;
        pred_xy_U3 = pred_xy_U1;
        
        pred_x_U6 = pred_x_U1;
        pred_y_U6 = pred_y_U1;
        pred_xy_U6 = pred_xy_U1;
        
        pred_x_U8 = pred_x_U1;
        pred_y_U8 = pred_y_U1;
        pred_xy_U8 = pred_xy_U1;
        
        %-----U2-U4-U5-U7 -----------
        pred_x_U2 = pred_x(SGs.ele_SGs/2+1:end,:,:);
        pred_y_U2 = pred_y(SGs.ele_SGs/2+1:end,:,:);
        pred_xy_U2 = pred_xy(SGs.ele_SGs/2+1:end,:,:);
        
        pred_x_U4 = pred_x_U2;
        pred_y_U4 = pred_y_U2;
        pred_xy_U4 = pred_xy_U2;
        
        pred_x_U5 = pred_x_U2;
        pred_y_U5 = pred_y_U2;
        pred_xy_U5 = pred_xy_U2;
        
        pred_x_U7 = pred_x_U2;
        pred_y_U7 = pred_y_U2;
        pred_xy_U7 = pred_xy_U2;
        
        if FE.dim ==2
            %---------------E12--------------
            pred_x_U6(:,3,:) = -pred_x_U6(:,3,:); pred_y_U6(:,3,:) = -pred_y_U6(:,3,:); pred_xy_U6(:,3,:) = -pred_xy_U6(:,3,:);
            pred_x_U8(:,3,:) = -pred_x_U8(:,3,:); pred_y_U8(:,3,:) = -pred_y_U8(:,3,:); pred_xy_U8(:,3,:) = -pred_xy_U8(:,3,:);
            
            pred_x_U5(:,3,:) = -pred_x_U5(:,3,:); pred_y_U5(:,3,:) = -pred_y_U5(:,3,:); pred_xy_U5(:,3,:) = -pred_xy_U5(:,3,:);
            pred_x_U7(:,3,:) = -pred_x_U7(:,3,:); pred_y_U7(:,3,:) = -pred_y_U7(:,3,:); pred_xy_U7(:,3,:) = -pred_xy_U7(:,3,:);
        else
            %---------------E12--------------
            pred_x_U6(:,4,:) = -pred_x_U6(:,4,:); pred_y_U6(:,4,:) = -pred_y_U6(:,4,:); pred_xy_U6(:,4,:) = -pred_xy_U6(:,4,:);
            pred_x_U8(:,4,:) = -pred_x_U8(:,4,:); pred_y_U8(:,4,:) = -pred_y_U8(:,4,:); pred_xy_U8(:,4,:) = -pred_xy_U8(:,4,:);
            
            pred_x_U5(:,4,:) = -pred_x_U5(:,4,:); pred_y_U5(:,4,:) = -pred_y_U5(:,4,:); pred_xy_U5(:,4,:) = -pred_xy_U5(:,4,:);
            pred_x_U7(:,4,:) = -pred_x_U7(:,4,:); pred_y_U7(:,4,:) = -pred_y_U7(:,4,:); pred_xy_U7(:,4,:) = -pred_xy_U7(:,4,:);
            %---------------E13--------------
            pred_x_U3(:,5,:) = -pred_x_U3(:,5,:); pred_y_U3(:,5,:) = -pred_y_U3(:,5,:); pred_xy_U3(:,5,:) = -pred_xy_U3(:,5,:);
            pred_x_U8(:,5,:) = -pred_x_U8(:,5,:); pred_y_U8(:,5,:) = -pred_y_U8(:,5,:); pred_xy_U8(:,5,:) = -pred_xy_U8(:,5,:);
            
            pred_x_U4(:,5,:) = -pred_x_U4(:,5,:); pred_y_U4(:,5,:) = -pred_y_U4(:,5,:); pred_xy_U4(:,5,:) = -pred_xy_U4(:,5,:);
            pred_x_U7(:,5,:) = -pred_x_U7(:,5,:); pred_y_U7(:,5,:) = -pred_y_U7(:,5,:); pred_xy_U7(:,5,:) = -pred_xy_U7(:,5,:);
            %---------------E23--------------
            pred_x_U3(:,6,:) = -pred_x_U3(:,6,:); pred_y_U3(:,6,:) = -pred_y_U3(:,6,:); pred_xy_U3(:,6,:) = -pred_xy_U3(:,6,:);
            pred_x_U6(:,6,:) = -pred_x_U6(:,6,:); pred_y_U6(:,6,:) = -pred_y_U6(:,6,:); pred_xy_U6(:,6,:) = -pred_xy_U6(:,6,:);
            
            pred_x_U4(:,6,:) = -pred_x_U4(:,6,:); pred_y_U4(:,6,:) = -pred_y_U4(:,6,:); pred_xy_U4(:,6,:) = -pred_xy_U4(:,6,:);
            pred_x_U5(:,6,:) = -pred_x_U5(:,6,:); pred_y_U5(:,6,:) = -pred_y_U5(:,6,:); pred_xy_U5(:,6,:) = -pred_xy_U5(:,6,:);            
            
        end
end
pred_x_Unit18 = [pred_x_U1;pred_x_U2;pred_x_U3;pred_x_U4;pred_x_U5;pred_x_U6;pred_x_U7;pred_x_U8];
pred_y_Unit18 = [pred_y_U1;pred_y_U2;pred_y_U3;pred_y_U4;pred_y_U5;pred_y_U6;pred_y_U7;pred_y_U8];
pred_xy_Unit18 = [pred_xy_U1;pred_xy_U2;pred_xy_U3;pred_xy_U4;pred_xy_U5;pred_xy_U6;pred_xy_U7;pred_xy_U8];

end
